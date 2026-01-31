using QuadGK
using StaticArrays

"""
    A_coeffs(b, c; rtol=1e-9) -> SVector{3,Float64}

Ellipsoid potential coefficients (starred): A* = a1 a2 a3 ∫ du / ((ai^2+u)Δ)
for (a1,a2,a3) = (1,b,c).
Sanity check: sum(A*) ≈ 2.
"""
function A_coeffs(b::Float64, c::Float64; rtol=1e-9)
    @assert 1.0 ≥ b ≥ c > 0.0
    integrand(u) = begin
        Δ = sqrt((1.0 + u) * (b^2 + u) * (c^2 + u))
        @SVector[
            1.0 / ((1.0 + u) * Δ),
            1.0 / ((b^2 + u) * Δ),
            1.0 / ((c^2 + u) * Δ),
        ]
    end
    val, _ = quadgk(integrand, 0.0, Inf; rtol=rtol)
    return (b * c) * val
end
