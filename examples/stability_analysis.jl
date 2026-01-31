include("../src/stability.jl")   # potential.jl も内部で include されます
include("../src/RocheEllipsoid.jl")  # ←ファイル名は適宜
using CairoMakie
using Makie.Colors
using .RocheEllipsoid
using Printf

deg(x) = 180 / pi * x

function fig3(ps)
    cmap = cgrad(:tab10; categorical=true)
    with_theme(theme_latexfonts()) do
        fig = Figure(fontsize=28)
        ax = Axis(fig[1, 1],
            xlabel=L"\arccos(a_3/a_1)\ \mathrm{[deg]}",
            ylabel=L"\sigma_3^2 / (\pi G \rho)",
        )
        open("Data/modes.dat", "w") do io
            @printf(io, "%8s%14s%8s%8s%8s%8s%8s%8s\n", "p", "acos(c)[deg]", "σ₁²", "σ₂²", "σ₃²", "σ₄²", "σ₅²", "σ₆²")
            for (i, p) in enumerate(ps)
                seq = roche_sequence(p; cmin=0.01, nc=500, rtolA=5e-9)   # (c,b,mu_over,omega_over)
                xs = Float64[]
                ys = Float64[]

                for (c, b, mu, om2) in zip(seq.c, seq.b, seq.mu_over, seq.omega_over)
                    σ_even = sort(even_sigma3sq(b, c, mu, om2; rtol=5e-9), rev=true)
                    σ_odd = sort(real.(odd_sigma2_roots(b, c, mu, om2; rtol=5e-9)), rev=true)
                    @printf(io, "%8.1e%14.2f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n",
                        p, deg(acos(c)),
                        σ_even[1], σ_even[2], σ_even[3],
                        σ_odd[1], σ_odd[2], σ_odd[3],
                    )

                    push!(xs, rad2deg(acos(c)))
                    push!(ys, minimum(σ_even))
                end
                # p>1e6の時はlabelをL"∞"にする
                # lbl = p >= 1e6 ? L"\infty" : "$(p)"
                lbl = @sprintf("%g", p)
                lines!(ax, xs, ys, label=lbl, color=cmap[i])
            end
        end
        axislegend(ax, L"p", position=:rt, framevisible=false)
        xlims!(ax, 0, 90)
        save("Plots/fig3_sigma3sq.png", fig)
        save("Plots/fig3_sigma3sq.pdf", fig)
        return fig
    end
end

ps = [-1.0, 0.0, 1.0, 4.0, 20.0, 100.0, 1e8]   # p→∞ は大きい数で近似
fig3(ps)
println("Saved: fig3_sigma3sq.png")
