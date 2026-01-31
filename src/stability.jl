# src/stability.jl
using LinearAlgebra
include("potential.jl")

# --- 小さなユーティリティ：3次方程式 roots（コンパニオン行列） ---
# s^3 + a2 s^2 + a1 s + a0 = 0
function cubic_roots_monic(a2, a1, a0)
    C = ComplexF64[
        -a2 -a1 -a0
        1 0 0
        0 1 0
    ]
    return eigvals(C)
end

# --- odd modes: 式(66)（λ^2 の3次）を解いて σ^2 = -λ^2 を返す ---
# λ^2(λ^2+4B13-Ω^2-μ)(λ^2+4B23-Ω^2+2μ) + 4Ω^2(λ^2+2B13+μ)(λ^2+2B23+μ)=0  :contentReference[oaicite:5]{index=5}
function odd_sigma2_roots(b::Float64, c::Float64, mu::Float64, Omega2::Float64; rtol=1e-9)
    B = B_star(b, c; rtol=rtol)
    Ω2 = Omega2

    α = 4B.B13 - Ω2 - mu
    β = 4B.B23 - Ω2 + 2mu
    γ = 2B.B13 + mu
    δ = 2B.B23 + mu

    # s = λ^2 として:
    # s(s+α)(s+β) + 4Ω^2(s+γ)(s+δ)=0
    # => s^3 + (α+β+4Ω^2)s^2 + (αβ + 4Ω^2(γ+δ))s + 4Ω^2 γδ = 0
    a2 = (α + β + 4Ω2)
    a1 = (α * β + 4Ω2 * (γ + δ))
    a0 = (4Ω2 * γ * δ)

    s_roots = cubic_roots_monic(a2, a1, a0) # λ^2 の根
    return -s_roots                         # σ^2 = -λ^2
end

# --- even modes: (67),(73),(74),(75) を cubic EVP にして λ を求める ---
# (67),(73),(74),(75) は page1195 の式 :contentReference[oaicite:6]{index=6}
function even_lambda_eigs(b::Float64, c::Float64, mu::Float64, Omega2::Float64; rtol=1e-9)
    B = B_star(b, c; rtol=rtol)
    ā = abar_axes(b, c)               # (75) の a_i は体積規格化を使う
    a1, a2, a3 = ā.a1, ā.a2, ā.a3

    Ω2 = Omega2
    Ω = sqrt(max(Ω2, 0.0))          # p=-1 で Ω=0 など

    α = 4B.B12 - 2Ω2 - mu            # (67) の定数部

    # unknown x = [V11, V22, V33, V12]
    n = 4
    M0 = zeros(ComplexF64, n, n)
    M1 = zeros(ComplexF64, n, n)
    M2 = zeros(ComplexF64, n, n)
    M3 = zeros(ComplexF64, n, n)

    # (67): (λ^2 + α)V12 + λΩ(V11 - V22)=0
    M2[1, 4] = 1
    M1[1, 1] = Ω
    M1[1, 2] = -Ω
    M0[1, 4] = α

    # (73):
    # (1/2 λ^2 - Ω^2 -2μ + 3B11 - B12) V11
    # - (1/2 λ^2 - Ω^2 + μ + 3B22 - B12) V22
    # + (B13 - B23) V33 - 2λΩ V12 =0
    M2[2, 1] = 1 / 2
    M2[2, 2] = -1 / 2
    M1[2, 4] = -2Ω
    M0[2, 1] = (-Ω2 - 2mu + 3B.B11 - B.B12)
    M0[2, 2] = (Ω2 - mu - 3B.B22 + B.B12)
    M0[2, 3] = (B.B13 - B.B23)

    # (74) は 1/λ があるので、式全体を λ 倍して cubic にする：
    # λE V11 + λF V22 + λG V33 + 6Ωμ V12 = 0
    # E = 1/2 λ^2 + (Ω^2 -2μ +3B11 +B12 -2B13)
    # F = 1/2 λ^2 + (Ω^2 + μ +3B22 +B12 -2B23)
    # G = -(λ^2 + (2μ +6B33 -B13 -B23))
    constE = (Ω2 - 2mu + 3B.B11 + B.B12 - 2B.B13)
    constF = (Ω2 + mu + 3B.B22 + B.B12 - 2B.B23)
    constG = (2mu + 6B.B33 - B.B13 - B.B23)

    M3[3, 1] = 1 / 2
    M3[3, 2] = 1 / 2
    M3[3, 3] = -1
    M1[3, 1] = constE
    M1[3, 2] = constF
    M1[3, 3] = -constG
    M0[3, 4] = 6Ω * mu

    # (75): V11/a1^2 + V22/a2^2 + V33/a3^2 = 0
    M0[4, 1] = 1 / a1^2
    M0[4, 2] = 1 / a2^2
    M0[4, 3] = 1 / a3^2

    # --- cubic polynomial eigenvalue: (λ^3 M3 + λ^2 M2 + λ M1 + M0)x=0 を線形化 ---
    Z = zeros(ComplexF64, n, n)
    Id = Matrix{ComplexF64}(LinearAlgebra.I, n, n)

    A = [Z Id Z
        Z Z Id
        -M0 -M1 -M2]

    Bmat = [Id Z Z
        Z Id Z
        Z Z M3]

    λs = eigvals(A, Bmat)
    # フィルタ（∞や巨大値を落とす）
    λs = [λ for λ in λs if isfinite(real(λ)) && isfinite(imag(λ)) && abs(λ) < 50]
    return λs
end

"""
even の σ^2 候補（σ^2 = -λ^2）から、
論文 Fig.3 の σ3^2 に相当する「最初に 0 を跨ぐ（=最小の）実数 σ^2」を返す。
"""
# 近い値を同一視して重複除去（ソート込み）
@inline function keep_canonical_lambda(λ::ComplexF64; tol=1e-12)
    re = real(λ)
    im = imag(λ)
    if im > tol
        return true
    elseif abs(im) ≤ tol && re > tol
        return true
    else
        return false
    end
end

function even_sigma3sq(b::Float64, c::Float64, mu::Float64, Omega2::Float64; rtol=1e-9)
    λs = even_lambda_eigs(b, c, mu, Omega2; rtol=rtol)
    isempty(λs) && return NaN

    cand = Float64[]
    for λ in λs
        keep_canonical_lambda(λ) || continue

        σ2 = -(λ^2)  # ← ここは括弧不要。σ2 = -(λ^2)

        re, im = real(σ2), imag(σ2)
        if isfinite(re) && isfinite(im) && abs(im) ≤ 1e-6 * (1 + abs(re))
            # 中立モード除外（ただし落としすぎ注意：必要なら eps→1e-12 などに）
            if abs(re) > eps()
                push!(cand, re)
            end
        end
    end

    isempty(cand) && return NaN
    return cand   # ← σ1^2 ≤ σ2^2 ≤ σ3^2 … の列
end
