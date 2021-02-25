include("general.jl")

## Types
# Impurity type
struct Impurity
    pos::Int    # Position of the impurity unit cell
    n::Int      # Index of the atom in the unit cell
    M::Float64  # Impurity mass in units of μ
    Δ::Float64  # External harmonic potential in units of k
end

# System type
struct System
    Ms::Vector{Float64}     # Masses in the 1D chain
    Imps::Vector{Impurity}  # Impurities in the system
    T::Float64              # Temperature of the system
    N::Int                  # Number of unit cells in the system (used for ED)
    K::Float64              # Confining potential for the atoms in the chain
end

## Exact Diagonalization 1D to obtain the total Free Energy
function exact_F(system)
    Ms = system.Ms
    Imps = system.Imps
    T = system.T
    N = system.N
    K = system.K
    # Number of atoms in the unit cell
    nAtom = length(Ms)
    if sum(map(x -> x.n, Imps) .<= nAtom) != length(Imps)
        error("Impurity index in the cell must not exceed the number of atoms")
    end
    # Prepare the matrix of masses
    M_List = repeat(Ms, N)
    # Prepare a list of unperturbed springs
    dv = (2 + K) .* ones(nAtom * N)
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
    Ω = sqrt.(abs.(Ω2))  # in units √(k / μ).
    # The free energy for each mode consists of the vacuum portion Ω / 2 and
    # the finite-T portion T * log(1 - exp(-Ω / T))
    total_energy = T * sum(log.(1 .- exp.(-Ω ./ T))) + sum(Ω) / 2
    return (total_energy)
end

# -TS term for the exact diagonalization
function exact_neg_TS(system)
    Ms = system.Ms
    Imps = system.Imps
    T = system.T
    N = system.N
    K = system.K
    # Number of atoms in the unit cell
    nAtom = length(Ms)
    if sum(map(x -> x.n, Imps) .<= nAtom) != length(Imps)
        error("Impurity index in the cell must not exceed the number of atoms")
    end

    # Prepare a list of unperturbed springs
    dv = (2 + K) .* ones(nAtom * N)
    # Prepare a list of Δs
    Δs = zeros(nAtom * N)

    for ii in Imps
        coord = nAtom * (ii.pos - 1) + ii.n
        Δs[coord] = ii.Δ
    end
    Δ_mat = Diagonal(Δs)
    # Assemble the matrices
    ev = -ones(nAtom * N - 1)
    U_Mat = SymTridiagonal(dv, ev) + zeros(nAtom * N, nAtom * N)
    # Make the system periodic by coupling the first and the last mass
    U_Mat[1, nAtom*N] = -1
    U_Mat[nAtom*N, 1] = -1

    res = T * log(det(Diagonal(ones(2 * N)) + inv(U_Mat) * Δ_mat)) / 2

    return res
end

function exact_F_I(system)
    Ms = system.Ms
    Imps = system.Imps
    T = system.T
    N = system.N
    K = system.K
    pristine_system = exact_F(System(Ms, [], T, N, K))
    total_energy = exact_F(system) - pristine_system
    impurity_energy =
        map(x -> exact_F(System(Ms, [x], T, N, K)) - pristine_system, Imps) |>
        sum
    return (total_energy - impurity_energy)
end

# -TS_I term for the exact diagonalization
function exact_neg_TS_I(system)
    Ms = system.Ms
    Imps = system.Imps
    T = system.T
    N = system.N
    K = system.K
    pristine_system = exact_zero_harmonic(System(Ms, [], T, N, K))
    total_energy = exact_zero_harmonic(system) - pristine_system
    impurity_energy =
        map(
            x ->
                exact_zero_harmonic(System(Ms, [x], T, N, K)) - pristine_system,
            Imps,
        ) |> sum
    return (total_energy - impurity_energy)
end

## Functions
@inline function modes(Ms, θ, K)
    # We set the spring constant k = 1
    # Inverse √[m] matrix
    M_Mat = Diagonal(1 ./ sqrt.(Ms))
    # Spring coupling matrix
    L_Mat = [
        2+K (-1-exp(-1im * θ))
        (-1-exp(1im * θ)) 2+K
    ]
    eig = eigen(Hermitian(M_Mat * L_Mat * M_Mat))
    return eig
end

@inline function Π_jl_integrand(z, Ms, K, θ, idx_j, idx_l)
    eig = modes(Ms, θ, K)
    freq = 1 ./ (-1 .+ eig.values ./ z^2)
    coupling = eig.vectors[idx_j, :] .* conj.(eig.vectors[idx_l, :])
    return (sum(coupling .* freq))
end

function Π_jl(z, Ms, K, Imp_j, Imp_l)
    D = Imp_j.pos - Imp_l.pos

    f_int(θ) = Π_jl_integrand(z, Ms, K, θ, Imp_j.n, Imp_l.n) * exp(1im * θ * D)

    res = quadgk(f_int, 0, 2 * π, atol = α, rtol = ν)[1]
    return (res / (2 * π * z^2))
end

function E_I_Integrand(ω, system)
    Ms = system.Ms
    imps = system.Imps
    K = system.K
    imp_host_idx = map(x -> x.n, imps)
    m = map(x -> Ms[x], imp_host_idx) |> Diagonal
    m_sqrt = sqrt(m)
    m_sqrt_inv = inv(m_sqrt)

    nImps = length(imps)
    imps_mat = repeat(imps, 1, nImps)
    imps_mat_T = permutedims(imps_mat)

    Π = map((x, y) -> Π_jl(ω, Ms, K, x, y), imps_mat, imps_mat_T)
    # Matrix to remove couplings between defects
    mask = Diagonal(ones(nImps))
    ΠD = mask .* Π

    Ξ = [
        m_sqrt_inv*Π*m_sqrt_inv ω*m_sqrt_inv*Π*m_sqrt
        ω*m_sqrt*Π*m_sqrt_inv (ω^2)*m_sqrt*Π*m_sqrt+m
    ]
    ΞD = [
        m_sqrt_inv*ΠD*m_sqrt_inv ω*m_sqrt_inv*ΠD*m_sqrt
        ω*m_sqrt*ΠD*m_sqrt_inv (ω^2)*m_sqrt*ΠD*m_sqrt+m
    ]

    unit_mat = ones(2 * length(imps)) |> Diagonal
    Δ_Λ =
        vcat(
            map(x -> x.Δ, imps),
            map(x -> 1 / x.M - 1 / system.Ms[x.n], imps),
        ) |> Diagonal
    return (
        (det(unit_mat .+ Ξ * Δ_Λ) |> Complex |> log) -
        (det(unit_mat .+ ΞD * Δ_Λ) |> Complex |> log)
    )
end


# E_I
function E_I(system)
    if system.T == 0
        res = quadgk(
            ω -> 2 * real(E_I_Integrand(ω * 1im, system)),
            0,
            Inf,
            maxevals = NumEvals,
            rtol = ν,
        )[1]::Float64
        return (res / (2 * π) / 2)
    else
        max_n = Integer(floor(max_omega / (2 * π * system.T)))
        res =
            map(
                n ->
                    system.T *
                    (E_I_Integrand(2 * pi * n * system.T * 1im, system)),
                1:1:max_n,
            ) |> sum |> real
        return res
    end
end

## Zero Matsubara Frequency: -TS
function Π0_jl_integrand(Ms, K, θ, idx_j, idx_l)
    eig = modes(Ms, θ, K)
    freq = 1 ./ eig.values
    coupling = eig.vectors[idx_j, :] .* conj.(eig.vectors[idx_l, :])
    return (sum(coupling .* freq))
end

function Π0_jl(Ms, K, Imp_j, Imp_l)
    D = Imp_j.pos - Imp_l.pos

    f_int(θ) = Π0_jl_integrand(Ms, K, θ, Imp_j.n, Imp_l.n) * exp(1im * θ * D)

    res = quadgk(f_int, 0, 2 * π, atol = α, rtol = ν)[1]
    return (res / (2 * π))
end

# -TS_I
function neg_TS_I(system)
    Ms = system.Ms
    imps = system.Imps
    K = system.K
    imp_host_idx = map(x -> x.n, imps)
    m = map(x -> Ms[x], imp_host_idx) |> Diagonal
    m_sqrt = sqrt(m)
    m_sqrt_inv = inv(m_sqrt)

    nImps = length(imps)
    imps_mat = repeat(imps, 1, nImps)
    imps_mat_T = permutedims(imps_mat)

    unit_mat = ones(length(imps)) |> Diagonal

    Π = map((x, y) -> Π0_jl(Ms, K, x, y), imps_mat, imps_mat_T)
    # Matrix to remove couplings between defects
    mask = Diagonal(ones(nImps))
    ΠD = mask .* Π

    res =
        (
            unit_mat +
            m_sqrt_inv * Π * m_sqrt_inv * Diagonal(map(x -> x.Δ, imps)) |>
            det |>
            log
        ) - (
            unit_mat +
            m_sqrt_inv * ΠD * m_sqrt_inv * Diagonal(map(x -> x.Δ, imps)) |>
            det |>
            log
        )
    return (res * system.T / 2) |> real
end
