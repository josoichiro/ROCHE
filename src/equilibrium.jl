using Roots

"""
    f_star(b, c, p; rtolA=1e-9)

(⋆) equation in cross-multiplied form (avoids 0/0 near the sphere limit).
"""
function f_star(b::Float64, c::Float64, p::Float64; rtolA=1e-9)
    A1, A2, A3 = A_coeffs(b, c; rtol=rtolA)  # A* assumed
    term1 = ((3.0 + p) + c^2) * (A2 * b^2 - A3 * c^2)
    term2 = (p * b^2 + c^2) * (A1 - A3 * c^2)
    return term1 - term2
end

"""
    solve_b_for_c(c, p; rtolA=1e-9, nscan=60) -> b or nothing

Bracket by scan on [c,1] and solve with bisection.
"""
function solve_b_for_c(c::Float64, p::Float64; rtolA=1e-9, nscan=60)
    bmin = max(c, 1e-12)
    bmax = 1.0
    f(b) = f_star(b, c, p; rtolA=rtolA)

    # endpoint root (sphere limit)
    fmax = try
        f(bmax)
    catch
        NaN
    end
    if isfinite(fmax) && abs(fmax) < 1e-12
        return bmax
    end

    bs = range(bmin, bmax; length=nscan)
    fs = [
        try
            f(b)
        catch
            NaN
        end for b in bs
    ]

    for i in 1:length(bs)-1
        fi, fj = fs[i], fs[i+1]
        if isfinite(fi) && isfinite(fj) && sign(fi) != sign(fj)
            return find_zero(f, (bs[i], bs[i+1]), Roots.Bisection())
        end
    end
    return nothing
end

"""
    roche_sequence(p; cmin=0.01, nc=160, rtolA=5e-9)

Scan c from 1 → cmin, solve b(c), and compute μ,Ω.
Returns NamedTuple with fields: c, b, mu_over, omega_over.
"""
function roche_sequence(p::Float64; cmin=0.01, nc=160, rtolA=5e-9)
    cs = range(1.0, cmin; length=nc)

    c_ok = Float64[]
    b_ok = Float64[]
    mu_over = Float64[]
    omega_over = Float64[]

    for c in cs
        b = solve_b_for_c(c, p; rtolA=rtolA)
        b === nothing && continue
        (1.0 ≥ b ≥ c) || continue

        μ, Ω = mu_omega_over(b, c, p; rtolA=rtolA)

        push!(c_ok, c)
        push!(b_ok, b)
        push!(mu_over, μ)
        push!(omega_over, Ω)
    end

    return (c=c_ok, b=b_ok, mu_over=mu_over, omega_over=omega_over)
end
