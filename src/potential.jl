# src/potential.jl
using QuadGK
using StaticArrays

"""
A* と Aij* を同時に返す（あなたのコードの A_coeffs は A* になっている前提）
(a1,a2,a3) = (1,b,c) を採用（比のみ必要なのでこれでOK）
戻り:
  Ai  :: SVector{3}  (A1*,A2*,A3*)
  Aij :: SMatrix{3,3} (Aij*; 対称)
"""
function A_Aij_star(b::Float64, c::Float64; rtol=1e-9)
    @assert 1.0 ≥ b ≥ c > 0.0
    a1, a2, a3 = 1.0, b, c
    pref = a1 * a2 * a3  # = b*c

    integrand(u) = begin
        Δ = sqrt((a1^2 + u) * (a2^2 + u) * (a3^2 + u))
        i1 = 1.0 / (a1^2 + u)
        i2 = 1.0 / (a2^2 + u)
        i3 = 1.0 / (a3^2 + u)

        # A_i* integrand: 1/((ai^2+u)Δ)
        A1 = i1 / Δ
        A2 = i2 / Δ
        A3 = i3 / Δ

        # A_ij* integrand: 1/((ai^2+u)(aj^2+u)Δ)
        A11 = (i1 * i1) / Δ
        A22 = (i2 * i2) / Δ
        A33 = (i3 * i3) / Δ
        A12 = (i1 * i2) / Δ
        A13 = (i1 * i3) / Δ
        A23 = (i2 * i3) / Δ

        @SVector [A1, A2, A3, A11, A22, A33, A12, A13, A23]
    end

    val, _ = quadgk(integrand, 0.0, Inf; rtol=rtol)
    v = pref * val

    Ai = @SVector [v[1], v[2], v[3]]
    Aij = @SMatrix [v[4] v[7] v[8];
        v[7] v[5] v[9];
        v[8] v[9] v[6]]
    return Ai, Aij
end

"""
B_ij* を返す（式(40): B_ij = A_i - a_j^2 A_ij）  
ここでは A*,Aij* を使うので同じ形でOK（B もスケール不変）。
戻り: NamedTuple (B11,B22,B33,B12,B13,B23)
"""
function B_star(b::Float64, c::Float64; rtol=1e-9)
    Ai, Aij = A_Aij_star(b, c; rtol=rtol)
    a = (1.0, b, c)
    # B(i,j) = Ai[i] - a[j]^2 * Aij[i,j]
    B11 = Ai[1] - a[1]^2 * Aij[1, 1]
    B22 = Ai[2] - a[2]^2 * Aij[2, 2]
    B33 = Ai[3] - a[3]^2 * Aij[3, 3]
    B12 = Ai[1] - a[2]^2 * Aij[1, 2]
    B13 = Ai[1] - a[3]^2 * Aij[1, 3]
    B23 = Ai[2] - a[3]^2 * Aij[2, 3]
    return (B11=B11, B22=B22, B33=B33, B12=B12, B13=B13, B23=B23)
end

"""
Fig.1/安定性の (75) に出てくる 1/a_i^2 用：体積規格化 ā_i（ā1ā2ā3=1）
ā1 = 1/(bc)^{1/3}, ā2=b/(bc)^{1/3}, ā3=c/(bc)^{1/3}  :contentReference[oaicite:4]{index=4}
"""
function abar_axes(b::Float64, c::Float64)
    s = (b * c)^(1 / 3)
    return (a1=1.0 / s, a2=b / s, a3=c / s)
end
