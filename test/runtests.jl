using HermiteUtils
using Test

@testset "HermiteUtils.jl" begin
    @test hermite_vandermonde(3) == 
        [1.0 1.0 1.0; 
         -1.2247448713915892 -8.881784197001252e-16 1.2247448713915892; 
         0.35355339059327395 -0.7071067811865475 0.35355339059327395];

    @test hermite_vandermonde(3, normalized=false) ==
        [1.0 1.0 1.0;
         -1.2247448713915892 -8.881784197001252e-16 1.2247448713915892;
         0.5000000000000003 -1.0 0.5000000000000003];

    @test He_n(3, 1.4) ≈ (1.4^3 - 3*1.4) / sqrt(6);
    @test He_n(3, -1.4, normalized=false) ≈ (-1.4^3 + 3*1.4);
    @test He_zero_to_n(3, 1.1) ≈ [1.0, 1.1, (1.1^2 - 1)/sqrt(2), (1.1^3 - 3.3)/sqrt(6)];
    @test He_zero_to_n(3, 1.1, normalized=false) ≈ [1.0, 1.1, 1.1^2 - 1, 1.1^3 - 3.3];
end
