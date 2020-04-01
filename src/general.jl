using LaTeXStrings
using LinearAlgebra
using Plots
using QuadGK

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

# System type
struct System
    Ms::Vector{Float64}     # Masses in the 1D chain
    Imps::Vector{Impurity}  # Impurities in the system
    T::Float64              # Temperature of the system
    N::Int                  # Number of unit cells in the system (used for ED)
end

## Exact Diagonalization 1D to obtain the total Free Energy
function exact_F(system)
    Ms = system.Ms
    Imps = system.Imps
    T = system.T
    N = system.N
    # Number of atoms in the unit cell
    nAtom = length(Ms)
    if sum(map(x -> x.n, Imps) .<= nAtom) != length(Imps)
        error("Impurity index in the cell must not exceed the number of atoms")
    end
    # Prepare a pristine chain potential energy matrix
    dv = 2 .* ones(nAtom * N)
    ev = -ones(nAtom * N - 1)
    U_Mat = SymTridiagonal(dv, ev) + zeros(nAtom * N, nAtom * N)
    # Make the system periodic by coupling the first and the last mass
    U_Mat[1, nAtom*N] = -1
    U_Mat[nAtom*N, 1] = -1
    # Prepare the matrix of masses
    M_List = repeat(Ms, N)
    # Replace the pristine masses by the impurities
    for ii in Imps
        coord = nAtom * (ii.pos - 1) + ii.n
        M_List[coord] = ii.λ
    end

    M_Mat = Diagonal(M_List)
    # Calculate the eigenvalues of the system which give
    # the squares of the energies
    Ω2 = eigvals(inv(M_Mat) * U_Mat) |> real
    Ω = sqrt.(abs.(Ω2[2:N*nAtom]))  # in units √(k / μ). We are dropping
    # the first mode because it has zero energy and can cause numerical issues

    # The free energy for each mode consists of the vacuum portion Ω / 2 and
    # the finite-T portion T * log(1 - exp(-Ω / T))
    total_energy = T * sum(log.(1 .- exp.(-Ω ./ T))) + sum(Ω) / 2
    return total_energy
end
## Functions
function modes(Ms, θ)
    M_Mat = Diagonal(1 ./ sqrt.(Ms))
    L_Mat = [
        2 (-1 - exp(-1im * θ))
        (-1 - exp(1im * θ)) 2
    ]
    eig = eigen(Hermitian(M_Mat * L_Mat * M_Mat))
    return eig
end

function Ξ_jl_integrand(z, Ms, θ, idx_j, idx_l)
    eig = modes(Ms, θ)
    freq = eig.values ./ (-z^2 .+ eig.values)
    coupling = eig.vectors[idx_j, :] .* conj.(eig.vectors[idx_l, :])
    return sum(coupling .* freq)
end

function Ξ_jl(z, system, j, l)
    Imp_j = system.Imps[j]
    Imp_l = system.Imps[l]

    mj = system.Ms[Imp_j.n]
    ml = system.Ms[Imp_l.n]
    D = Imp_j.pos - Imp_l.pos |> abs

    f_int(θ) =
        -√(mj * ml) *
        Ξ_jl_integrand(z, system.Ms, θ, Imp_j.n, Imp_l.n) *
        exp(1im * θ * D)

    res = quadgk(f_int, 0, 2 * π)[1]
    return (res / (2 * π))
end

function Ξ(z, system)
    nImps = length(system.Imps)
    js = repeat(1:nImps, 1, nImps)
    ls = transpose(js)
    res = map((j, l) -> Ξ_jl(z, system, j, l), js, ls)
    return res
end

# Interaction energy (for now computed at T = 0)
# FI Integrand for zero-T energy
function F_I_Integrand_T0(system, z)
    nImps = length(system.Imps)
    Ξ_ = Ξ(z, system)
    λs = map(x -> x.λ, system.Imps)
    ms = map(x -> system.Ms[x.n], system.Imps)
    α = Diagonal(1 ./ ms .- 1 ./ λs)
    Δ = (Matrix{Int}(I, nImps, nImps) .+ Ξ_ * α)
    Ξ0 = map(j -> Ξ_jl(z, system, j, j), 1:nImps) |> Diagonal
    Δ0_Inv = (Matrix{Int}(I, nImps, nImps) .+ Ξ0 * α) |> inv

    return (det(Δ0_Inv * Δ) |> Complex |> log |> real)
end

# FI for zero T
function F_I_T0(system)
    f(w) = F_I_Integrand_T0(system, 1im * w)
    res = quadgk(f, -Inf, 0, Inf, maxevals = NumEvals)[1]
    return (res / (2 * π))
end

## TODO
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
my_red = RGB(215 / 255, 67 / 255, 84 / 255)
my_green = RGB(106 / 255, 178 / 255, 71 / 255)
my_blue = RGB(100 / 255, 101 / 255, 218 / 255)
my_violet = RGB(169 / 255, 89 / 255, 201 / 255)
my_orange = RGB(209 / 255, 135 / 255, 46 / 255)

colors = [my_red, my_green, my_blue, my_violet, my_orange]
