module RocheEllipsoid

export A_coeffs,
    f_star,
    solve_b_for_c,
    mu_omega_over,
    roche_sequence,
    abar12_from_seq

include("./coeffs.jl")
include("./equilibrium.jl")
include("./derived.jl")

end # module
