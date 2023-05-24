# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Dates
using TaylorIntegration
using Test

using InteractiveUtils: methodswith

@testset "Integration methods" begin
    
    @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_ng_eph!}, TaylorIntegration.jetcoeffs!))
    @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_ng_eph!}, TaylorIntegration._allocate_jetcoeffs!))

    @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_ng_eph_threads!}, TaylorIntegration.jetcoeffs!))
    @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_ng_eph_threads!}, TaylorIntegration._allocate_jetcoeffs!))

    @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_eph_threads!}, TaylorIntegration.jetcoeffs!))
    @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_eph_threads!}, TaylorIntegration._allocate_jetcoeffs!))
    
end 