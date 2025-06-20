using NEOs
using Dates
using Test

@testset "Common" begin

    @testset "Time conversions" begin
        t0 = now()
        @test et2dtutc(dtutc2et(t0)) == t0
        @test jdtdb2dtutc(dtutc2jdtdb(t0)) == t0
    end

end