include("../src/RocheEllipsoid.jl")
using CairoMakie
using Makie.Colors
using Printf
using .RocheEllipsoid

function main()
    ps = [-1.0, 0.0, 1.0, 4.0, 20.0, 100.0, 1e8]
    seqs = Dict(p => roche_sequence(p; cmin=0.01, nc=500, rtolA=5e-9) for p in ps)
    # ---- Fig.1-like: (ā1, ā2) ----
    cmap = cgrad(:tab10; categorical=true)
    with_theme(theme_latexfonts()) do
        fig1 = Figure(fontsize=28)
        ax1 = Axis(fig1[1, 1],
            xlabel=L"\bar{a}_1",
            ylabel=L"\bar{a}_2",
        )

        for (i, p) in enumerate(ps)
            lbl = @sprintf("%g", p)
            s = seqs[p]
            isempty(s.c) && continue
            a1, a2 = abar12_from_seq(s)
            idx = sortperm(a1)
            lines!(ax1, a1[idx], a2[idx], label=lbl, color=cmap[i])
        end
        axislegend(ax1, L"p", position=:rt, framevisible=false)
        xlims!(ax1, 1.0, 2.3)
        ylims!(ax1, 0.6, 1.3)
        save("Plots/roche_fig1_like.png", fig1)
        save("Plots/roche_fig1_like.pdf", fig1)
    end

    # ---- Fig.2: Ω²/(πGρ) vs φ ----
    with_theme(theme_latexfonts()) do
        fig2 = Figure(fontsize=28)
        ax2 = Axis(fig2[1, 1],
            xlabel=L"\arccos(a_3/a_1)\ \mathrm{[deg]}",
            ylabel=L"\Omega^2/(\pi G \rho)",
        )
        for (i, p) in enumerate(ps)
            lbl = @sprintf("%g", p)
            s = seqs[p]
            isempty(s.c) && continue
            φ = rad2deg.(acos.(s.c))
            lines!(ax2, φ, s.omega_over, label=lbl, color=cmap[i])
            imax = argmax(s.omega_over)
            scatter!(ax2, [φ[imax]], [s.omega_over[imax]], color=cmap[i])
            @printf("p=%g: max at φ=%.2f deg, Ω²/(πGρ)=%.4f\n", p, φ[imax], s.omega_over[imax])
        end
        axislegend(ax2, L"p", position=:lt, framevisible=false)
        xlims!(ax2, 0.0, 90.0)
        save("Plots/roche_fig2_like.png", fig2)
        save("Plots/roche_fig2_like.pdf", fig2)
    end

    # println("Saved: roche_fig1_like.png, roche_fig2_like.png")
end

main()
