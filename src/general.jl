using Plots
using QuadGK
using LinearAlgebra

## Parameters
const ν = 1e-3;         # Relative tolerance for integration
const η = 1e-7;         # Small number used for i0
const NumEvals = 1e7;   # Maximum number of evaluations in quadgk

# Impurity type
struct Impurity
    pos::Int    # Position of the impurity unit cell
    n::Int      # Index of the atom in the unit cell
    λ::Float64  # Impurity mass in units of μ
end

## Functions
function Q2_sθ(Ms, s, θ)
    m = Ms[1]
    M = Ms[2]
    # s should be 1 for optical branch, 0 for acoustic branch
    # - eventually generalise to n branches
    sqTerm = sqrt(m^2 + M^2 + 2 * m * M * cos(θ)) / (m * M)
    invariantTerm = (m + M) / (m * M)
    if s == 0
        return (invariantTerm - sqTerm)
    elseif s == 1
        return (invariantTerm + sqTerm)
    else
        error("This branch is not available at the moment")
    end
end

# Unnormalized mode vector
function e_sθ(Ms, s, θ)
    m = Ms[1]
    M = Ms[2]
    return ([1 + exp(-1im * θ), 2 - m * Q2_sθ(m, M, s, θ)])
end

# Mode couling term
function coupling(Ms, s, θ)
    e = e_sθ(Ms, s, θ)
    return (Diagonal(Ms) * e * e' * Diagonal(Ms)) / (e' * Diagonal(Ms) * e)
end

# Impurity coupling
function YGY(Ms, z, D, nj, nk)
    f_int_acc(θ) =
        Q2_sθ(Ms, 0, θ) * exp(1im * D * θ) / (-z^2 + Q2_sθ(Ms, 0, θ)) *
        (coupling(Ms, 0, θ)[nj, nk])

    f_int_opt(θ) =
        Q2_sθ(Ms, 1, θ) * exp(1im * D * θ) / (-z^2 + Q2_sθ(Ms, 1, θ)) *
        (coupling(Ms, 1, θ)[nj, nk])

    f_int(θ) = f_int_opt(θ) + f_int_acc(θ)
    res = quadgk(f_int, 0, 2 * π)[1]
    return (-res) / (2 * π)
end

function Δ(Ms, z, Imps)
    nImps = length(Imps)
    α = map(x -> 1 / Ms[x.n] - 1 / x.λ, Imps) |> Diagonal
    Imp_Mat = repeat(Imps, 1, nImps)    # Impurity position matrix
    ImpT_Mat = permutedims(Imp_Mat)     # Transpose of impurity position matrix
    YGY_ = map(
        (x, y) -> YGY(m, M, z, abs.(x.pos - y.pos), x.n, y.n),
        Imp_Mat,
        ImpT_Mat,
    )

    return (Matrix{Int}(I, nImps, nImps) .+ Prop_Mat * α)
end

function FI_Integrand(Ms, z, Imps)
    α = map(x -> 1 / Ms[x.n] - 1 / x.λ, Imps) |> Diagonal
    Δ_ = Δ(Ms, z, Imps)
    Δ0_Inv =
        map(
            x -> 1 ./ (1 .+ YGY(Ms, z, 0, x.n, x.n) .* (1 / Ms[x.n] - 1 / x.λ)),
            Imps,
        ) |> Diagonal

    return det(Δ0_Inv * Δ_) |> Complex |> log |> real
end

@time YGY([1, 2], 1.2 + 1im * η, 2, 1,2)

# Q2_sθ([1,2.], 1, 1.2)
