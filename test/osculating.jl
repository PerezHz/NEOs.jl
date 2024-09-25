using NEOs
using Dates
using Test
using LinearAlgebra: norm

@testset "Keplerian osculating elements" begin
    q = [-9.759018085743707E-01, 3.896554445697074E-01, 1.478066121706831E-01, -9.071450085084557E-03, -9.353197026254517E-03, -5.610023032269034E-03]
    kep = pv2kep(q)
    jd2000 = datetime2julian(DateTime(2000,1,1,12))
    @test norm(kep(jd2000) - q, Inf) < 1e-10

    qa = [-1.0506628055913627, -0.06064314196134998, -0.04997102228887035, 0.0029591421121582077, -0.01423233538611057, -0.005218412537773594]
    jda = datetime2julian(DateTime(2004,6,1))
    kepa = pv2kep(qa)
    @test norm(kepa(jd2000) - qa, Inf) < 1e-10
end