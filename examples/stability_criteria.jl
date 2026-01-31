using DelimitedFiles
using Printf
using Roots
using Interpolations

data = readdlm("Data/modes.dat", skipstart=1)
ps = unique(data[:, 1])
@printf("%6s %8s\n", "p", "φDI[deg]")
for p in ps
    idxs = findall(data[:, 1] .== p)
    acos_c = data[idxs, 2]
    sigma3sq = data[idxs, 5]

    itp = LinearInterpolation(acos_c, sigma3sq, extrapolation_bc=Line())
    degs = 30.0:0.1:70.0
    vals = itp(degs)
    # root finding for σ₃² = 0
    f(d) = itp(d)
    d_root = find_zero(f, (30.0, 70.0))
    @printf("%6.0e %8.2f\n", p, d_root)
    # for (d, v) in zip(degs, vals)
    #     @printf("%8.2f%14.6f\n", d, v)
    # end
end