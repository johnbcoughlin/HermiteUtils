module HermiteUtils

using FastGaussQuadrature
using LinearAlgebra

export hermitepoints, hermiteweights, hermite_vandermonde, weighted_hermite_vandermonde,
    He_n, He_zero_to_n, hermiteanalysis, hermitesynthesis

"""
    hermitepoints(n)

Returns a vector of the Gauss-Hermite quadrature points of order n.
These are exactly the roots of He_n.
"""
hermitepoints(n) = gausshermite(n)[1] .* sqrt(2)

hermiteweights(n) = gausshermite(n)[2] / sqrt(π)

"""
    hermite_vandermonde(n, [x = hermitepoints(n)])

Returns the Vandermonde matrix whose entries are He_{i-1}(x_j) / (i-1)!,
in the ith row and jth column.

The default value of x is hermitepoints(n).
"""
function hermite_vandermonde(n; normalized=true)
    hermite_vandermonde(n, hermitepoints(n); normalized)
end

function hermite_vandermonde(n, x; normalized=true)
    V = Array{Float64}(undef, length(x), n)
    sqrt_factorial = vcat(1.0, cumprod(1 .* sqrt.(1:n-1)))
    for k in eachindex(x)
        V[k, :] = FastGaussQuadrature.hermpoly_rec(0:n-1, x[k]) * exp(x[k]^2/4)
        if !normalized
            V[k, :] .*= sqrt_factorial
        end
    end
    V'
end

function weighted_hermite_vandermonde(n; normalized=true)
    V = Array{Float64}(undef, n, n)
    x = hermitepoints(n)

    sqrt_factorial = vcat(1.0, cumprod(1 .* sqrt.(1:n-1)))

    for k in 1:n
        V[k, :] = FastGaussQuadrature.hermpoly_rec(0:n-1, sqrt(2) * x[k])
        if !normalized
            V[k, :] .*= sqrt_factorial
        end
    end
    V'
end

"""
    hermiteanalysis(n)

Returns a matrix which computes the coefficients f_i in the following expansion:

f = ∑ exp(-x^2/2) He_i(x) f_i

The function values f must be evaluated at hermitepoints(n).
"""
function hermiteanalysis(n)
    x = hermitepoints(n)
    w = hermiteweights(n)

    hermite_vandermonde(n) * Diagonal(w .* exp.(x.^2 ./ 2))
end

hermitesynthesis(n) = hermitesynthesis(n, hermitepoints(n))

function hermitesynthesis(n, x)
    Diagonal(exp.(-x.^2 ./ 2)) * hermite_vandermonde(n, x)'
end

"""
    He_n(n, x; normalized=true)

Computes the vector [He_0(x), He_1(x), ..., He_n(x)]
"""
function He_zero_to_n(n, x; normalized=true)
    sqrt_factorial = vcat(1.0, cumprod(1 .* sqrt.(1:n)))
    result = FastGaussQuadrature.hermpoly_rec(0:n, x) * exp(x^2/4)
    if normalized
        result
    else
        result .* sqrt_factorial
    end
end

"""
    He_n(n, x; normalized=true)

Computes He_n(x)
"""
function He_n(n::Int, x; normalized=true)
    sqrt_factorial = vcat(1.0, cumprod(1 .* sqrt.(1:n)))
    result = FastGaussQuadrature.hermpoly_rec(0:n, x)[end] * exp(x^2/4)
    if normalized
        result
    else
        result .* sqrt_factorial[end]
    end
end


end
