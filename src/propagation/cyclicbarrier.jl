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

# Split a block into its statements, each paired with the LineNumberNode that
# precedes it (or `nothing`), so the generated code keeps the source lines
function lines_of(ex::Expr)
    lines = Tuple{Any, Any}[]
    lnn = nothing
    for x in ex.args
        if x isa LineNumberNode
            lnn = x
        else
            push!(lines, (lnn, x))
            lnn = nothing
        end
    end
    return lines
end

# Flatten a vector of (LineNumberNode, statement) pairs
linenums(lnn) = isnothing(lnn) ? Any[] : Any[lnn]
unlines(lines) = mapreduce(l -> Any[linenums(l[1])..., l[2]], vcat, lines; init = Any[])

# Remove the LineNumberNodes that do not point to `file`, e.g. those coming
# from the quotes in this file or from the expansion of Threads.@threads
filter_linenums!(x, file) = x
function filter_linenums!(ex::Expr, file)
    if ex.head === :block || ex.head === :quote
        filter!(x -> !(x isa LineNumberNode && x.file !== file), ex.args)
    end
    foreach(x -> filter_linenums!(x, file), ex.args)
    return ex
end

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

# Assemble the (LineNumberNode, statement) pairs into blocks
function assemble_blocks(x::AbstractVector)
    mask = [is_threads_macro(l[2]) for l in x]
    iters = Vector{Any}(undef, 0)
    blocks = Vector{Vector{Tuple{Any, Any}}}(undef, 0)
    for i in eachindex(x)
        if mask[i]
            push!(blocks, [x[i]])
            push!(iters, x[i][2].args[end].args[1].args[2])
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
    if is_threads_macro(x[1][2])
        lnn, ex = x[1]
        loop = ex.args[end]
        loop_def, loop_body = loop.args
        itr = Symbol(:chunked_indices_, i)
        chunk = Expr(:for, :($(loop_def.args[1]) = $(itr)[task_id]), loop_body)
        q = Expr(:block, linenums(lnn)..., chunk, :(wait($cyclic_barrier)))
    else
        serial = Expr(:if, :(task_id == 1), Expr(:block, unlines(x)...))
        q = Expr(:block, serial, :(wait($cyclic_barrier)))
    end
    return q
end

# Build the multi-threaded (cyclic barrier) version of a vector of
# (LineNumberNode, statement) pairs (unescaped). Also return the symbols
# assigned in the serial code outside the main loop
function cyclicbarrier_threaded(lines::AbstractVector)
    # Split the lines into the code before, within and after the main loop
    kmain = findall(l -> is_main_loop(l[2]), lines)
    length(kmain) <= 1 || throw(ArgumentError("At most one top-level for loop \
        can contain multi-threaded loops"))
    if isempty(kmain)
        prelines, mainline, postlines = lines, nothing, lines[1:0]
    else
        k = kmain[1]
        prelines, mainline, postlines = lines[1:k-1], lines[k], lines[k+1:end]
    end
    # Assemble the lines of code into blocks
    preiters, preblocks = assemble_blocks(prelines)
    if isnothing(mainline)
        main_lnn, loop_def = nothing, nothing
        bodyiters, bodyblocks = assemble_blocks(lines[1:0])
    else
        main_lnn, mainloop = mainline
        loop_def, loop_body = mainloop.args
        bodyiters, bodyblocks = assemble_blocks(lines_of(loop_body))
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
    main = isnothing(loop_def) ? nothing :
        Expr(:block, linenums(main_lnn)..., Expr(:for, loop_def, Expr(:block, body...)))
    # Variables assigned in the serial code outside the main loop must be
    # shared by all tasks and visible after the block
    shared = Set{Symbol}()
    for b in vcat(preblocks, postblocks)
        is_threads_macro(b[1][2]) && continue
        foreach(l -> assigned_symbols!(shared, l[2]), b)
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
    return sort!(collect(shared)), q
end

# Fields of `DynamicalParameters` holding an `EphemerisEvaluationBuffer`
const EPHEMERIS_BUFFERS = (:sseph, :acceph, :poteph)

# Check if an expression is `local x = params.eph(t)`, where `eph` is one of
# the `EPHEMERIS_BUFFERS`
function is_ephemeris_evaluation(ex)
    ex isa Expr && ex.head === :local && length(ex.args) == 1 || return false
    asg = ex.args[1]
    asg isa Expr && asg.head === :(=) && asg.args[1] isa Symbol || return false
    call = asg.args[2]
    call isa Expr && call.head === :call && length(call.args) == 2 || return false
    f = call.args[1]
    return f isa Expr && f.head === :. && f.args[1] isa Symbol &&
           f.args[2] isa QuoteNode && f.args[2].value in EPHEMERIS_BUFFERS
end

# Unfold `local x = params.eph(t)` into (i) the `local` declarations needed to
# evaluate the ephemeris and (ii) the `Threads.@threads` evaluation loop
function unfold_ephemeris(ex::Expr)
    x, call = ex.args[1].args
    p, eph = call.args[1].args[1], call.args[1].args[2].value
    tt = call.args[2]
    t_, e_, a_, T_, i_, δ_ = Symbol.(eph, (:_t, :_eph, :_aux, :_ephT, :_ind, :_δt))
    decls = Any[
        :(local $t_ = $p.$eph.t),
        :(local $e_ = $p.$eph.eph),
        :(local $a_ = $p.$eph.aux),
        :(local $T_ = $p.$eph.ephT),
        :(local $x = $p.$eph.ephU),
        :(TaylorSeries.identity!($t_, $tt, 0)),
        Expr(:local, Expr(:(=), Expr(:tuple, i_, δ_), :(timeindex($e_, $t_)))),
    ]
    loop = quote
        Threads.@threads for i in eachindex($x)
            TaylorSeries.zero!($T_[i])
            TaylorSeries.zero!($a_[i])
            TaylorSeries._horner!($T_[i], $e_.p[$i_, i], $δ_, $a_[i])
            for k in eachindex($x[i])
                taylorembed!($x[i], $T_[i], k)
            end
        end
    end
    Base.remove_linenums!(loop)
    return decls, loop.args[1]
end

# Kind of function decorated by @cyclicbarrier
function cyclicbarrier_kind(fname)
    fname == :(TaylorIntegration.jetcoeffs!) && return :jetcoeffs
    fname == :(TaylorIntegration._allocate_jetcoeffs!) && return :allocate
    fname isa Symbol && return :model
    throw(ArgumentError("@cyclicbarrier cannot decorate function $fname"))
end

# Name of the serial / multi-threaded methods of each kind of function
cyclicbarrier_inner(kind::Symbol, fname) = kind === :jetcoeffs ? :_jetcoeffs! :
    kind === :allocate ? :_allocate! : fname

# Give a name to every argument of a function signature
function name_arguments(args)
    names, newargs = Symbol[], Any[]
    for (i, a) in enumerate(args)
        if a isa Symbol
            push!(names, a)
            push!(newargs, a)
        elseif a isa Expr && a.head === :(::) && length(a.args) == 2 && a.args[1] isa Symbol
            push!(names, a.args[1])
            push!(newargs, a)
        elseif a isa Expr && a.head === :(::) && length(a.args) == 1
            name = Symbol(:__arg, i)
            push!(names, name)
            push!(newargs, Expr(:(::), name, a.args[1]))
        else
            throw(ArgumentError("@cyclicbarrier does not support argument $a"))
        end
    end
    return names, newargs
end

# Build a function signature
function build_signature(fname, args, wparams)
    call = Expr(:call, fname, args...)
    return isempty(wparams) ? call : Expr(:where, call, wparams...)
end

is_local(x) = x isa Expr && x.head === :local

# Generate the wrapper, serial and multi-threaded methods of a function (unescaped)
function cyclicbarrier_function(fdef::Expr, source::LineNumberNode)
    # Split the function definition
    sig, body = fdef.args
    wparams = Any[]
    if sig isa Expr && sig.head === :where
        append!(wparams, sig.args[2:end])
        sig = sig.args[1]
    end
    @assert sig isa Expr && sig.head === :call "Unsupported function signature"
    fname, args = sig.args[1], sig.args[2:end]
    kind = cyclicbarrier_kind(fname)
    inner = cyclicbarrier_inner(kind, fname)
    names, newargs = name_arguments(args)
    :params in names || throw(ArgumentError("@cyclicbarrier requires an argument \
        named `params`"))
    # Split the function's body into: (i) the header, i.e. the lines before the
    # first `local` declaration, (ii) the `local` declarations, (iii) the rest
    # of the code and (iv) the final `return` statement (if any)
    lines = lines_of(body)
    tail = !isempty(lines) && lines[end][2] isa Expr && lines[end][2].head === :return ?
        [pop!(lines)] : lines[1:0]
    k = findfirst(l -> is_local(l[2]), lines)
    k = isnothing(k) ? length(lines) + 1 : k
    header, rest = lines[1:k-1], lines[k:end]
    locals = filter(l -> is_local(l[2]), rest)
    code = filter(l -> !is_local(l[2]), rest)
    # Unfold the evaluation of the Solar System ephemerides
    decls, ephloops = lines[1:0], lines[1:0]
    for (lnn, line) in locals
        if is_ephemeris_evaluation(line)
            d, loop = unfold_ephemeris(line)
            append!(decls, [(lnn, x) for x in d])
            push!(ephloops, (lnn, loop))
        else
            push!(decls, (lnn, line))
        end
    end
    # Serial version
    unthreaded(l) = (l[1], remove_threads(l[2]))
    serial = Expr(:block, unlines(header)..., unlines(decls)...,
        unlines(unthreaded.(ephloops))..., unlines(unthreaded.(code))..., unlines(tail)...)
    # Multi-threaded version
    if kind === :model
        # Non-parsed dynamical models keep their Threads.@threads loops
        threaded = Expr(:block, unlines(header)..., unlines(decls)...,
            unlines(ephloops)..., unlines(code)..., unlines(tail)...)
    else
        # Parsed methods use a cyclic barrier
        shared, q = cyclicbarrier_threaded(vcat(ephloops, code))
        shared_decls = [Expr(:local, s) for s in shared]
        threaded = Expr(:block, unlines(header)..., unlines(decls)..., shared_decls...,
            Expr(:let, Expr(:block), q), unlines(tail)...)
    end
    # Wrapper with the original signature
    wsig = build_signature(fname, newargs, wparams)
    wbody = quote
        if params.threads && Threads.threadpoolsize() > 1
            return $inner(Val(true), $(names...))
        else
            return $inner(Val(false), $(names...))
        end
    end
    pushfirst!(wbody.args, source)
    ssig = build_signature(inner, Any[:(::Val{false}), newargs...], wparams)
    tsig = build_signature(inner, Any[:(::Val{true}), newargs...], wparams)
    q = quote
        Base.@__doc__ $(Expr(:function, wsig, wbody))
        $(Expr(:function, ssig, serial))
        $(Expr(:function, tsig, threaded))
    end
    return q
end

"""
    @cyclicbarrier function f(...) ... end

Generate serial and multi-threaded versions of a dynamical model function
(`f(dq, q, params, t)`) or of its parsed methods of
`TaylorIntegration._allocate_jetcoeffs!` and `TaylorIntegration.jetcoeffs!`
generated by `@taylorize`. Three methods are defined:

- A wrapper, with the original signature, that calls the multi-threaded version
    if `params.threads && Threads.threadpoolsize() > 1`, and the serial version
    otherwise.
- The serial version, where every `Threads.@threads` loop is run as a plain
    `for` loop.
- The multi-threaded version. For the parsed methods, the `Threads.@threads`
    loops are run via a cyclic barrier, so tasks are spawned only once per call.
    Dynamical model functions keep their `Threads.@threads` loops.

The serial and multi-threaded versions are named `_jetcoeffs!` / `_allocate!`
for the parsed methods, and `f` for dynamical models, with an additional first
argument `::Val{false}` / `::Val{true}`. The function's body must be
organized as follows:

- An optional header, i.e. the lines before the first `local` declaration
    (e.g. `order = ...` and the unpacking of `__ralloc` in `jetcoeffs!`),
- `local` declarations, which are hoisted before the rest of the code,
- the rest of the code, which (for `jetcoeffs!`) may contain at most one
    top-level `for` loop with `Threads.@threads` loops in its body, and
- an optional final `return` statement.

Every `local x = params.eph(t)` declaration, where `eph` is one of `:sseph`,
`:acceph` or `:poteph`, is unfolded into the explicit evaluation of the
corresponding `EphemerisEvaluationBuffer`, whose (multi-threaded) evaluation
loop is placed at the beginning of the code.

To keep the generated code readable, its nested macros are expanded and only the
`LineNumberNode`s pointing to the decorated function's file are kept; these are
needed by code coverage tools and stack traces.

!!! warning
    This macro is on an experimental stage; check the integration results carefully.
"""
macro cyclicbarrier(fdef)
    (fdef isa Expr && fdef.head === :function) || throw(ArgumentError("@cyclicbarrier \
        must decorate a function definition"))
    q = cyclicbarrier_function(fdef, __source__)
    # Expand the nested macros (e.g. Threads.@threads and Threads.@spawn) and keep
    # only the LineNumberNodes of the decorated function's file, so the generated
    # code is easier to read while code coverage and stack traces still work
    q = macroexpand(__module__, q; recursive = true)
    filter_linenums!(q, __source__.file)
    # We use esc() to ensure variables resolve in the caller's scope (macro hygiene)
    return esc(q)
end
