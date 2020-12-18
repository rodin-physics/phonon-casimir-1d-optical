using LaTeXStrings
using LinearAlgebra
using Plots
using QuadGK

## Parameters
const ν = 1e-5;         # Relative tolerance for integration
const η = 1e-5;         # Small number used for i0
const α = 1e-9;         # Small number for absolute tolerance
const max_omega = 10000; # large number used to truncate the finite temperature sum
const NumEvals = 1e7;   # Maximum number of evaluations in quadgk

# Impurity type
struct Impurity
    pos::Int    # Position of the impurity unit cell
    n::Int      # Index of the atom in the unit cell
    M::Float64  # Impurity mass in units of μ
    Δ::Float64  # External harmonic potential in units of K
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
    # Prepare the matrix of masses
    M_List = repeat(Ms, N)
    # Prepare a list of unperturbed springs
    dv = 2 .* ones(nAtom * N)
    # Replace the pristine masses by the impurities
    # and add a harmonic confinement term
    for ii in Imps
        coord = nAtom * (ii.pos - 1) + ii.n
        M_List[coord] = ii.M
        dv[coord] = dv[coord] + ii.Δ
    end

    # Assemble the matrices
    ev = -ones(nAtom * N - 1)
    U_Mat = SymTridiagonal(dv, ev) + zeros(nAtom * N, nAtom * N)
    # Make the system periodic by coupling the first and the last mass
    U_Mat[1, nAtom*N] = -1
    U_Mat[nAtom*N, 1] = -1
    M_Mat = Diagonal(M_List)
    # Calculate the eigenvalues of the system which give
    # the squares of the energies
    Ω2 = eigvals(inv(M_Mat) * U_Mat) |> real
    Ω = sqrt.(abs.(Ω2[2:end]))  # in units √(k / μ). We are dropping
    # the first mode because it has zero energy and can cause numerical issues

    # The free energy for each mode consists of the vacuum portion Ω / 2 and
    # the finite-T portion T * log(1 - exp(-Ω / T))
    total_energy = T * sum(log.(1 .- exp.(-Ω ./ T))) + sum(Ω) / 2
    return total_energy
end
## Functions
@inline function modes(Ms, θ)
    # We set the spring constant K = 1
    # Inverse √[m] matrix
    M_Mat = Diagonal(1 ./ sqrt.(Ms))
    # Spring coupling matrix
    L_Mat = [
        2 (-1-exp(-1im * θ))
        (-1-exp(1im * θ)) 2
    ]
    eig = eigen(Hermitian(M_Mat * L_Mat * M_Mat))
    return eig
end

function Π_jl_integrand(z, Ms, θ, idx_j, idx_l)
    eig = modes(Ms, θ)
    freq = 1 ./ (-1 .+ eig.values ./ z^2)
    coupling = eig.vectors[idx_j, :] .* conj.(eig.vectors[idx_l, :])
    return (sum(coupling .* freq))
end

function Π_jl(z, Ms, Imp_j, Imp_l)
    D = Imp_j.pos - Imp_l.pos

    f_int(θ) = Π_jl_integrand(z, Ms, θ, Imp_j.n, Imp_l.n) * exp(1im * θ * D)

    res = quadgk(f_int, 0, 2 * π, atol = α, rtol = ν)[1]
    return (res / (2 * π * z^2))
end

function Ξ(z, system)
    Ms = system.Ms
    imps = system.Imps
    imp_host_idx = map(x -> x.n, imps)
    m = map(x -> Ms[x], imp_host_idx) |> Diagonal
    m_sqrt = sqrt(m)
    m_sqrt_inv = inv(m_sqrt)

    nImps = length(imps)
    imps_mat = repeat(imps, 1, nImps)
    imps_mat_T = permutedims(imps_mat)

    Π = map((x, y) -> Π_jl(z, Ms, x, y), imps_mat, imps_mat_T)

    res = [
        m_sqrt_inv*Π*m_sqrt_inv z*m_sqrt_inv*Π*m_sqrt
        z*m_sqrt*Π*m_sqrt_inv (z^2)*m_sqrt*Π*m_sqrt+m
    ]
    return (res)

end
#
function F_I_Integrand_T0(ω, system)
    Ms = system.Ms
    imps = system.Imps
    imp_host_idx = map(x -> x.n, imps)
    m = map(x -> Ms[x], imp_host_idx) |> Diagonal
    m_sqrt = sqrt(m)
    m_sqrt_inv = inv(m_sqrt)

    nImps = length(imps)
    imps_mat = repeat(imps, 1, nImps)
    imps_mat_T = permutedims(imps_mat)

    Π = map((x, y) -> Π_jl(1im * ω, Ms, x, y), imps_mat, imps_mat_T)

    Ξ = [
        m_sqrt_inv*Π*m_sqrt_inv (1im*ω)*m_sqrt_inv*Π*m_sqrt
        (1im*ω)*m_sqrt*Π*m_sqrt_inv ((1im*ω)^2)*m_sqrt*Π*m_sqrt+m
    ]

    unit_mat = ones(2 * length(imps)) |> Diagonal
    Δ_Λ =
        vcat(map(x -> x.Δ, imps), map(x -> 1 / x.M - 1 / system.Ms[x.n], imps)) |> Diagonal
    return (det(unit_mat .+ Ξ * Δ_Λ) |> Complex |> log |> real)
end

function F_I_Integrand(ω, system)
    T = system.T
    Ms = system.Ms
    imps = system.Imps
    imp_host_idx = map(x -> x.n, imps)
    m = map(x -> Ms[x], imp_host_idx) |> Diagonal
    m_sqrt = sqrt(m)
    m_sqrt_inv = inv(m_sqrt)

    nImps = length(imps)
    imps_mat = repeat(imps, 1, nImps)
    imps_mat_T = permutedims(imps_mat)

    Π = map((x, y) -> Π_jl(ω, Ms, x, y), imps_mat, imps_mat_T)

    Ξ = [
        m_sqrt_inv*Π*m_sqrt_inv ω*m_sqrt_inv*Π*m_sqrt
        ω*m_sqrt*Π*m_sqrt_inv (ω^2)*m_sqrt*Π*m_sqrt+m
    ]

    unit_mat = ones(2 * length(imps)) |> Diagonal
    Δ_Λ =
        vcat(map(x -> x.Δ, imps), map(x -> 1 / x.M - 1 / system.Ms[x.n], imps)) |> Diagonal
    return imag(log(det(unit_mat + Ξ * Δ_Λ))) / (exp(ω / T) - 1)
end


function F_I_Integrand_Discrete(ω, system)
    T = system.T
    Ms = system.Ms
    imps = system.Imps
    imp_host_idx = map(x -> x.n, imps)
    m = map(x -> Ms[x], imp_host_idx) |> Diagonal
    m_sqrt = sqrt(m)
    m_sqrt_inv = inv(m_sqrt)

    nImps = length(imps)
    imps_mat = repeat(imps, 1, nImps)
    imps_mat_T = permutedims(imps_mat)

    Π = map((x, y) -> Π_jl(ω, Ms, x, y), imps_mat, imps_mat_T)

    Ξ = [
        m_sqrt_inv*Π*m_sqrt_inv ω*m_sqrt_inv*Π*m_sqrt
        ω*m_sqrt*Π*m_sqrt_inv (ω^2)*m_sqrt*Π*m_sqrt+m
    ]

    unit_mat = ones(2 * length(imps)) |> Diagonal
    Δ_Λ =
        vcat(map(x -> x.Δ, imps), map(x -> 1 / x.M - 1 / system.Ms[x.n], imps)) |> Diagonal
    return log(det(unit_mat + Ξ * Δ_Λ))
end

function F_I_T0(system)
    res =
        quadgk(
            ω -> 2 * real(F_I_Integrand_T0(ω, system)),
            0,
            max_omega,
            maxevals = NumEvals,
            rtol = ν,
        )[1]::Float64
    return (res / (2 * π) / 2) # See the F_I formula for the division by 2
end


function F_I_T_Discrete(system)
    max_n = Integer(floor(max_omega / (2 * π * system.T)))
    res =
        map(
            n -> system.T * F_I_Integrand_Discrete(2 * pi * n * system.T * 1im, system),
            range(1, max_n, length = max_n),
        ) |> sum
    return (res) # See the F_I formula for the division by 2
end

function F_I(system)
    res =
        quadgk(
            x -> real(F_I_Integrand(x + 1im * η, system)),
            -Inf,
            0,
            Inf,
            maxevals = NumEvals,
            rtol = ν,
        )[1]::Float64
    return res / (2 * pi)
end


# Colors for plotting
my_red = RGB(215 / 255, 67 / 255, 84 / 255)
my_green = RGB(106 / 255, 178 / 255, 71 / 255)
my_blue = RGB(100 / 255, 101 / 255, 218 / 255)
my_violet = RGB(169 / 255, 89 / 255, 201 / 255)
my_orange = RGB(209 / 255, 135 / 255, 46 / 255)

colors = [my_red, my_green, my_blue, my_violet, my_orange]
