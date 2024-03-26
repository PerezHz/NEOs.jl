using Test
using NEOs
using Aqua

@testset "Aqua tests (performance)" begin
    # This tests that we don't accidentally run into
    # https://github.com/JuliaLang/julia/issues/29393
    # Aqua.test_unbound_args(NEOs)
    ua = Aqua.detect_unbound_args_recursively(NEOs)
    @test length(ua) == 0

    # See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/1750
    # Test that we're not introducing method ambiguities across deps
    ambs = Aqua.detect_ambiguities(NEOs; recursive = true)
    pkg_match(pkgname, pkdir::Nothing) = false
    pkg_match(pkgname, pkdir::AbstractString) = occursin(pkgname, pkdir)
    filter!(x -> pkg_match("NEOs", pkgdir(last(x).module)), ambs)
    for method_ambiguity in ambs
        @show method_ambiguity
    end
    if VERSION < v"1.10.0-DEV"
        @test length(ambs) == 0
    end
end

@testset "Aqua tests (additional)" begin
    Aqua.test_undefined_exports(NEOs)
    Aqua.test_deps_compat(NEOs)
    Aqua.test_stale_deps(NEOs)
    Aqua.test_piracies(NEOs)
    Aqua.test_unbound_args(NEOs)
    Aqua.test_project_extras(NEOs)
    Aqua.test_persistent_tasks(NEOs)
end

