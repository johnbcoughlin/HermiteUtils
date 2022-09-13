using HermiteUtils
using Test
using LinearAlgebra

@testset "HermiteUtils.jl" begin
    @test norm(He_n.(10, hermitepoints(10))) < sqrt(eps())

    @test hermite_vandermonde(3) == 
        [0.9999999999999999 1.0 0.9999999999999999; 
         -1.7320508075688774 -1.2560739669470201e-15 1.7320508075688774; 
         1.4142135623730956 -0.7071067811865475 1.4142135623730956];

    @test hermite_vandermonde(3, normalized=false) ==
        [0.9999999999999999 1.0 0.9999999999999999; 
         -1.7320508075688774 -1.2560739669470201e-15 1.7320508075688774; 
         2.000000000000001 -1.0 2.000000000000001]

    @test He_n(3, 1.4) ≈ (1.4^3 - 3*1.4) / sqrt(6);
    @test He_n(3, -1.4, normalized=false) ≈ (-1.4^3 + 3*1.4);
    @test He_zero_to_n(3, 1.1) ≈ [1.0, 1.1, (1.1^2 - 1)/sqrt(2), (1.1^3 - 3.3)/sqrt(6)];
    @test He_zero_to_n(3, 1.1, normalized=false) ≈ [1.0, 1.1, 1.1^2 - 1, 1.1^3 - 3.3];


    x = hermitepoints(100)
    e4 = zeros(100);
    e4[4] = 1.0
    ha = hermiteanalysis(100)
    hs = hermitesynthesis(100)

    fM = exp.(-x.^2 ./ 2)
    f = (x.^3 .- 3x) .* fM

    @test ha * f ≈ sqrt(6) * e4
    @test hs * sqrt(6) * e4 ≈ f
end
