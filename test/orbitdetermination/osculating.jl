using NEOs
using Dates
using Test
using LinearAlgebra: norm

@testset "Keplerian osculating elements" begin
    q = [-9.759018085743707E-01, 3.896554445697074E-01, 1.478066121706831E-01,
         -9.071450085084557E-03, -9.353197026254517E-03, -5.610023032269034E-03]
    kep = pv2kep(q)
    jd = datetime2julian(DateTime(2000, 1, 1, 12))
    @test norm(kep(jd) - q, Inf) < 1e-10

    q = [-1.0506628055913627, -0.06064314196134998, -0.04997102228887035,
          0.0029591421121582077, -0.01423233538611057, -0.005218412537773594]
    kep = pv2kep(q)
    @test norm(kep(jd) - q, Inf) < 1e-10
end