using Plots
using QuadGK
using LinearAlgebra
using LaTeXStrings

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

## Exact Diagonalization 1D
function exact_F(Ms, N, Imps, T)
    # Number of atoms in the unit cell
    nAtom = length(Ms)
    # Prepare a pristine chain potential energy matrix
    dv = 2 .* ones(nAtom * N)
    ev = -ones(nAtom * N - 1)
    U_Mat = SymTridiagonal(dv, ev) + zeros(nAtom * N, nAtom * N)
    # Make the system periodic by coupling the first and the last mass
    U_Mat[1, nAtom*N] = -1
    U_Mat[nAtom*N, 1] = -1
    # Prepare the matrix of masses
    M_Mat = repeat(Ms, N)
    # Replace the pristine masses by the impurities
    for ii in Imps
        coord = nAtom * (ii.pos - 1) + ii.n
        M_Mat[coord] = ii.λ
    end
    M_Mat = Diagonal(M_Mat)
    # Calculate the eigenvalues of the system which give
    # the squares of the energies
    Ω2 = eigvals(inv(M_Mat) * U_Mat)
    Ω = sqrt.(abs.(Ω2[2:N*nAtom]))  # in units √(k / μ). We are dropping
    # the first mode because it has zero energy

    # The free energy for each mode consists of the vacuum portion Ω / 2 and
    # the finite-T portion T * log(1 - exp(-Ω / T))
    total_energy = T * sum(log.(1 .- exp.(-Ω ./ T))) + sum(Ω) / 2
    return total_energy
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
    return ([1 + exp(-1im * θ), 2 - m * Q2_sθ([m, M], s, θ)])
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

# @time YGY([1, 2], 1.2 + 1im * η, 2, 1, 2)
# YGY([1, 2], 1.2 + 1im * η, 2, 1, 2)
# # Q2_sθ([1,2.], 1, 1.2)
#
#
# exact_F([1], 40, [Impurity(1, 1, 2.0), Impurity(4, 1, 2.0)], 1e-2)
# exact_F([1, 1], 20, [Impurity(1, 1, 2.0), Impurity(2, 2, 2.0)], 1e-2)
# Exact_Free_Energy(5, 2.0, 2.0, 1, 2 * 1e-2) * 2


# Colors for plotting
my_red = RGB(215/255,67/255,84/255)
my_green = RGB(106/255,178/255,71/255)
my_blue = RGB(100/255,101/255,218/255)
my_violet = RGB(169/255,89/255,201/255)
my_orange = RGB(209/255,135/255,46/255)

colors = [my_red
        , my_green
        , my_blue
        , my_violet
        , my_orange
         ]
