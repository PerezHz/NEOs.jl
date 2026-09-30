"""
    CyclicBarrier

Synchronization tool used to make a set of tasks wait for each
other to reach a common barrier point in the code.

# Fields

- `n::Int`: number of tasks to wait for.
- `count::Int`: number of tasks that have reached the barrier.
- `cond::Threads.Condition`: trigger to wake up all waiting tasks.
- `generation::Int`: number of times the barrier has been used.
"""
mutable struct CyclicBarrier
    n::Int
    count::Int
    cond::Threads.Condition
    generation::Int
end

# Outer constructor
CyclicBarrier(n::Int) = CyclicBarrier(n, 0, Threads.Condition(), 0)

# Override Base.wait
function wait(x::CyclicBarrier)
    lock(x.cond) do
        gen = x.generation
        x.count += 1
        # The last task to arrive resets the count, increments the generation,
        # and wakes up all waiting tasks
        if x.count == x.n
            x.count = 0
            x.generation += 1
            notify(x.cond)
        # All other tasks wait until the generation number changes
        else
            while gen == x.generation
                wait(x.cond)
            end
        end
    end
end

# Replicate the distribution of a loop's iterations over threads
# done by Julia's built-in static scheduler
function static_loop_split(r::AbstractVector{Int}, N::Int, tid::Int)
    len = length(r)
    chunk_size = div(len, N)
    remainder = rem(len, N)
    start_offset = (tid - 1) * chunk_size + min(tid - 1, remainder)
    num_elements = chunk_size + (tid <= remainder ? 1 : 0)
    return (first(r) + start_offset):(first(r) + start_offset + num_elements - 1)
end

# Check if an expression is a call to Threads.@threads
is_threads_macro(x) = false
function is_threads_macro(ex::Expr)
    # Check if the expression is a macro call
    ex.head === :macrocall || return false
    # The first argument of a :macrocall is the macro name
    macro_name = ex.args[1]
    # Match against common ways @threads is invoked
    return macro_name === Symbol("@threads") ||
           macro_name == :(Threads.@threads) ||
           macro_name == :(Threads.var"@threads") ||
           macro_name == :(Base.Threads.@threads)
end

# Replace every Threads.@threads loop in an expression by a plain for loop
remove_threads(x) = x
function remove_threads(ex::Expr)
    is_threads_macro(ex) && return remove_threads(ex.args[end])
    return Expr(ex.head, map(remove_threads, ex.args)...)
end

# Filter out LineNumberNodes
clean_blocks(x::Expr) = map(identity, filter(!Base.Fix2(isa, LineNumberNode), x.args))

# Check if an expression is the main loop, i.e. a for loop with at least
# one Threads.@threads loop at the top level of its body
is_main_loop(x) = x isa Expr && x.head === :for && any(is_threads_macro, clean_blocks(x.args[2]))

# Collect the symbols (re)bound by assignments within an expression
assigned_symbols!(s::Set{Symbol}, x) = s
function assigned_symbols!(s::Set{Symbol}, ex::Expr)
    # Do not descend into nested functions (they have their own scope)
    ex.head in (:function, :->) && return s
    # Skip the iteration specification of for loops
    ex.head === :for && return assigned_symbols!(s, ex.args[2])
    if ex.head in (:(=), :+=, :-=, :*=, :/=, :^=)
        lhs = ex.args[1]
        # Short-form function definitions
        lhs isa Expr && lhs.head in (:call, :where) && return s
        lhs_symbols!(s, lhs)
        return assigned_symbols!(s, ex.args[2])
    end
    foreach(Base.Fix1(assigned_symbols!, s), ex.args)
    return s
end

lhs_symbols!(s::Set{Symbol}, x) = s
lhs_symbols!(s::Set{Symbol}, x::Symbol) = push!(s, x)
function lhs_symbols!(s::Set{Symbol}, x::Expr)
    x.head === :tuple && foreach(Base.Fix1(lhs_symbols!, s), x.args)
    return s
end

# Assemble the lines of code into blocks
function assemble_blocks(x::AbstractVector)
    mask = is_threads_macro.(x)
    iters = Vector{Any}(undef, 0)
    blocks = Vector{Vector{Any}}(undef, 0)
    for i in eachindex(x)
        if mask[i]
            push!(blocks, [x[i]])
            push!(iters, x[i].args[end].args[1].args[2])
        else
            if isempty(blocks) || mask[i-1]
                push!(blocks, [x[i]])
                push!(iters, nothing)
            else
                push!(blocks[end], x[i])
            end
        end
    end
    return iters, blocks
end

function modify_block(x::AbstractVector, i::Int)
    cyclic_barrier = Symbol(:cyclic_barrier_, i)
    if is_threads_macro(x[1])
        loop = x[1].args[end]
        loop_def, loop_body = loop.args
        itr = Symbol(:chunked_indices_, i)
        q = quote
            for $(loop_def.args[1]) in $(itr)[task_id]
                $(loop_body.args...)
            end
            wait($cyclic_barrier)
        end
    else
        q = quote
           if task_id == 1
               $(x...)
           end
           wait($cyclic_barrier)
        end
    end
    Base.remove_linenums!(q)
    return q
end

# Build the multi-threaded (cyclic barrier) version of a block (unescaped). Also
# return the symbols assigned in the serial code outside the main loop
function cyclicbarrier_threaded(ex::Expr)
    # Filter out LineNumberNodes to find the actual lines of code
    lines = clean_blocks(ex)
    # Split the lines into the code before, within and after the main loop
    kmain = findall(is_main_loop, lines)
    length(kmain) <= 1 || throw(ArgumentError("At most one top-level for loop \
        can contain multi-threaded loops"))
    if isempty(kmain)
        prelines, mainloop, postlines = lines, nothing, lines[1:0]
    else
        k = kmain[1]
        prelines, mainloop, postlines = lines[1:k-1], lines[k], lines[k+1:end]
    end
    # Assemble the lines of code into blocks
    preiters, preblocks = assemble_blocks(prelines)
    if isnothing(mainloop)
        loop_def = nothing
        bodyiters, bodyblocks = assemble_blocks(lines[1:0])
    else
        loop_def, loop_body = mainloop.args
        bodyiters, bodyblocks = assemble_blocks(clean_blocks(loop_body))
    end
    postiters, postblocks = assemble_blocks(postlines)
    iters = vcat(preiters, bodyiters, postiters)
    blocks = vcat(preblocks, bodyblocks, postblocks)
    mblocks = modify_block.(blocks, eachindex(blocks))
    # Preamble, body and postamble
    npre, nbody = length(preblocks), length(bodyblocks)
    getargs(x) = mapreduce(Base.Fix2(getfield, :args), vcat, x; init = [])
    preamble = getargs(view(mblocks, 1:npre))
    body = getargs(view(mblocks, npre+1:npre+nbody))
    postamble = getargs(view(mblocks, npre+nbody+1:length(mblocks)))
    # Reconstruct the original main loop
    main = isnothing(loop_def) ? nothing : quote
        for $(loop_def.args[1]) in $(loop_def.args[2])
            $(body...)
        end
    end
    # Variables assigned in the serial code outside the main loop must be
    # shared by all tasks and visible after the block
    shared = Set{Symbol}()
    for b in vcat(preblocks, postblocks)
        is_threads_macro(b[1]) && continue
        foreach(Base.Fix1(assigned_symbols!, shared), b)
    end
    # Variable declarations
    names = Symbol.(:cyclic_barrier_, eachindex(mblocks))
    cyclic_barriers = [:($name = CyclicBarrier(Ntasks)) for name in names]
    names = Symbol.(:chunked_indices_, eachindex(iters))
    chunked_indices = [:($name = [static_loop_split($iter, Ntasks, tid) for tid in 1:Ntasks])
        for (name, iter) in zip(names, iters) if !isnothing(iter)]
    # Generate the multithreaded code
    q = quote
        # Variable declarations
        Ntasks = Threads.threadpoolsize()
        $(cyclic_barriers...)
        $(chunked_indices...)
        function threadsfor_fun(task_id::Int)
            # Preamble
            $(preamble...)
            # Main loop
            $main
            # Postamble
            $(postamble...)
        end
        # Spawn tasks
        tasks = Vector{Task}(undef, Ntasks)
        for task_id in 1:Ntasks
            tasks[task_id] = Threads.@spawn threadsfor_fun(task_id)
        end
        # Wait for all tasks to complete
        for t in tasks
            fetch(t)
        end
    end
    Base.remove_linenums!(q)
    return sort!(collect(shared)), q
end

function cyclicbarrier_expr(threads, ex)
    # Verify the provided expression is actually a block
    @assert ex isa Expr && ex.head === :block "Expression must be a block"
    shared, threaded = cyclicbarrier_threaded(ex)
    serial = remove_threads(ex)
    decls = [Expr(:local, s) for s in shared]
    # We use esc() to ensure variables resolve in the caller's scope (macro hygiene)
    q = esc(quote
        $(decls...)
        if $threads
            let
                $threaded
            end
        else
            $serial
        end
    end)
    Base.remove_linenums!(q)
    return q
end

"""
    @cyclicbarrier [threads] ex

This macro modifies a block of code in the `jetcoeffs!` and `_allocate_jetcoeffs!`
functions generated by `@taylorize` in order to implement a cyclic barrier that
avoids spawning new tasks at every `Threads.@threads` loop. Tasks are spawned
only once per evaluation of the block, and each `Threads.@threads` loop is
replaced by a statically-scheduled chunk of iterations followed by a barrier.

The block may contain:
- `Threads.@threads` loops (e.g. the Solar System ephemeris evaluation loops),
- serial code, which is run by the first task only, and
- at most one top-level `for` loop (the main loop, e.g. the loop over the orders
    of the Taylor expansions in `jetcoeffs!`) whose body contains `Threads.@threads`
    loops at its top level.

Variables assigned in the serial code outside the main loop are declared `local`
in the enclosing scope, so they are shared by all tasks and remain visible after
the block. The iteration ranges of the `Threads.@threads` loops must be computable
before the block is entered.

If `threads` (a `Bool` expression evaluated at runtime, e.g. `params.threads`)
is `false`, every `Threads.@threads` loop inside `ex` is run as a plain
serial `for` loop instead. If omitted, `threads` defaults to `true`.

!!! warning
    This macro is on an experimental stage; check the integration results carefully.
"""
macro cyclicbarrier(threads, ex)
    return cyclicbarrier_expr(threads, ex)
end

macro cyclicbarrier(ex)
    return cyclicbarrier_expr(true, ex)
end

"""
    @optionalthreads threads [schedule] for ... end

Run a `for` loop multi-threaded via `Threads.@threads` (with the optional
`schedule` argument, e.g. `:static`) if `threads` is `true`, and as a plain
serial `for` loop otherwise. `threads` is a `Bool` expression evaluated at
runtime, e.g. `params.threads && Threads.threadpoolsize() > 1`.

This macro is used in the non-parsed dynamical models of
`src/propagation/dynamicalmodels.jl`. Since `@taylorize` only accepts
`Threads.@threads`, every `@optionalthreads threads` must be replaced by
`Threads.@threads` before generating the parsed methods of `jetcoeffs!`
(see the header of `src/propagation/jetcoeffs.jl`).
"""
macro optionalthreads(threads, args...)
    isempty(args) && throw(ArgumentError("@optionalthreads must be followed \
        by a `for` loop"))
    loop = args[end]
    (loop isa Expr && loop.head === :for) || throw(ArgumentError("@optionalthreads \
        must be followed by a `for` loop"))
    threaded = Expr(:macrocall, Expr(:., :Threads, QuoteNode(Symbol("@threads"))),
        __source__, args...)
    # We use esc() to ensure variables resolve in the caller's scope (macro hygiene)
    return esc(quote
        if $threads
            $threaded
        else
            $loop
        end
    end)
end
