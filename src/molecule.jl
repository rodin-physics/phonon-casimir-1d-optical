include("general.jl")
## Exact Diagonalization for the Diatomic Molecule
function Ω_molecule(δ1, δ2, m1, m2, K)
    mass_mat = [m1 0; 0 m2]
    force_mat = [1+K+δ1 (-1); (-1) 1+K+δ2]
    Ω_sq = eigvals(inv(mass_mat) * force_mat)
    if K == 0
        Ω_sq = maximum(Ω_sq)
    end
    return sqrt.(Ω_sq)
end

function F_exact_molecule(δ1, δ2, m1, m2, K, T)
    Ωs = Ω_molecule(δ1, δ2, m1, m2, K)
    total_energy = T * sum(log.(1 .- exp.(-Ωs ./ T))) + sum(Ωs) / 2
    return total_energy
end

function F_I_exact_molecule(δ, M, K, T)
    res =
        F_exact_molecule(δ, δ, M, M, K, T) -
        2 * F_exact_molecule(δ, 0, M, 1, K, T) +
        F_exact_molecule(0, 0, 1, 1, K, T)
    return res
end

## QFT Energy
function E_I_Integrand_molecule(ω, δ, M, K, T)
    Ω1_sq = K
    Ω2_sq = 2 + K
    Π_diag = (1 / (-ω^2 + Ω1_sq) + 1 / (-ω^2 + Ω2_sq)) / 2
    Π_off_diag = (1 / (-ω^2 + Ω1_sq) - 1 / (-ω^2 + Ω2_sq)) / 2
    Π = [Π_diag Π_off_diag; Π_off_diag Π_diag]
    ΠD = [Π_diag 0; 0 Π_diag]

    mh = [1 0; 0 1]
    zrs = zeros(2, 2)
    left_vec = [inv(mh) .^ (1 / 2); ω * mh .^ (1 / 2)]
    Ξ = left_vec * Π * transpose(left_vec) + [zrs zrs; zrs mh]
    ΞD = left_vec * ΠD * transpose(left_vec) + [zrs zrs; zrs mh]

    Δ = [δ 0; 0 δ]
    Λ = [(1/M-1) 0; 0 (1/M-1)]
    Δ_Λ = [Δ zrs; zrs Λ]

    one_mat = Diagonal(ones(4))
    res = log(det(one_mat + Ξ * Δ_Λ)) - log(det(one_mat + ΞD * Δ_Λ))
    return res

end

function E_I_molecule(δ, M, K, T)
    if T == 0
        res =
            quadgk(x -> E_I_Integrand_molecule(1im * x, δ, M, K, T), 0, Inf)[1] /
            (2 * π)

    else
        res = sum(
            map(
                n ->
                    T * E_I_Integrand_molecule(2im * π * T * n, δ, M, K, T),
                1:10000,
            ),
        )
    end
    return real.(res)
end

function neg_TS_I_molecule(δ, M, K, T)
    Ω1_sq = K
    Ω2_sq = 2 + K
    Π_diag = (1 / (Ω1_sq) + 1 / (Ω2_sq)) / 2
    Π_off_diag = (1 / (Ω1_sq) - 1 / (Ω2_sq)) / 2
    Π = [Π_diag Π_off_diag; Π_off_diag Π_diag]
    ΠD = [Π_diag 0; 0 Π_diag]
    one_mat = Diagonal(ones(2))
    mh = sqrt.(inv([1 0; 0 1]))
    Δ = [δ 0; 0 δ]
    res =
        log(det(one_mat + mh * Π * mh * Δ)) -
        log(det(one_mat + mh * ΠD * mh * Δ))
    return real.(T * res / 2)
end

function F_I_molecule(δ, M, K, T)
    return (E_I_molecule(δ, M, K, T) + neg_TS_I_molecule(δ, M, K, T))
end


# F_I_molecule(2, 1, 1e-1, 0)
