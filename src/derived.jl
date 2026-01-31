"""
    mu_omega_over(b, c, p; rtolA=1e-9) -> (mu/(πGρ), Omega^2/(πGρ))

Uses starred A* returned by A_coeffs.
"""
function mu_omega_over(b::Float64, c::Float64, p::Float64; rtolA=1e-9)
    A1, A2, A3 = A_coeffs(b, c; rtol=rtolA)
    mu_over = 2.0 * (A1 - A3 * c^2) / ((3.0 + p) + c^2)
    omega_over = (1.0 + p) * mu_over
    return mu_over, omega_over
end

"""
    abar12_from_seq(seq) -> (ā1, ā2)

Volume-normalized axes for a1=1, a2=b, a3=c:
ā1 = 1/(bc)^{1/3}, ā2 = b/(bc)^{1/3}.
"""
function abar12_from_seq(seq)
    b = seq.b
    c = seq.c
    s = (b .* c) .^ (1 / 3)
    abar1 = 1.0 ./ s
    abar2 = b ./ s
    return abar1, abar2
end
