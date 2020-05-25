using QuadGK
using Plots
using LinearAlgebra
pyplot()

const errtol = 1e-9;
T = 1e-12;
const NumEvals = 1e7;
const ν = 1e-4;         # Relative tolerance for integration
const η = 1e-7;         # Small number used for i0

# Impurity type
struct ImpuritySL
    p::Int      # Index of the unit cell in p
    λ::Float64  # Impurity mass in units of μ
end

# System type
struct SystemSL
    M::Float64    # Mass of an atom
    Imps::Vector{ImpuritySL}  # Impurities in the system
    k::Float64              # longitudinal spring constant
    α::Float64              # transversal spring constant
    T::Float64              # Temperature of the system
    N::Int                  # Length of system (sqrt number of unit cells)
end

# first three functions for navigating and placing stuff in the lattice
function p2nm(N, p)
    if p > N ^ 2
        p = p % N ^ 2
    end
    n = Int(floor((p - 1) / N)) + 1
    m = (p - ( n - 1 ) * N) % N
    if m == 0
        m = N
    end

    return n, m
end

function nm2p(N, n, m)

    n = n % N
    m = m % N
    if n == 0
        n = N
    end
    if m == 0
        m = N
    end

    p = (((n - 1) * N) + m ) % (N ^ 2)

    if p == 0
        p = N ^ 2
    end

    return p
end

function get_adj(N, p, movealong, pm1)
    # movealong can be 'r' for row and 'c' for column
    # pm1 is plus or minus 1 for whether we moving in positive or negative direction

    n, m = p2nm(N, p)
    if movealong == 'r'
        # move along the row, ie left and right
        newp = nm2p(N, n, m + pm1)
    elseif movealong == 'c'
        # move along the column, ie up and down
        newp = nm2p(N, n + pm1, m)
    end

    return newp
end


function A2(θ, k)
    return k^2 + sin(θ / 2)^2
end

function asoln1(θ, ω, Ω, D2)
    k = ω / Ω
    Asq = A2(θ, k)
    b = 2 * Asq + 1
    b2m1 = b ^ 2 - 1
    return 4 * pi * (b - sqrt(b2m1)) ^ D2 / sqrt(b2m1)
end

function asoln2(θ, ω, Ω, D2)
    k = ω / Ω
    Asq = A2(θ, k)
    b = 2 * Asq + 1
    b2m1 = b ^ 2 - 1
    rm = b - sqrt(b2m1)

    factor0 = pi
    factor1 = rm ^ D2
    factor2 = (2 - rm - 1 / rm)
    factor3 = 1 / sqrt(b2m1)
    return factor0 * factor1 * factor2 * factor3
end

function integral1(ω, Ω, D1, D2)
    f_int1(θ) = (1 / (4 * pi ^ 2)) * exp( 1im * D1 * θ) * sin(θ / 2) * asoln1(θ, ω, Ω, D2)
    res = quadgk(f_int1, 0, 2 * pi, atol = errtol)[1]
    return res
end

function integral2(ω, Ω, D1, D2)
    f_int2(θ) = (1 / (4 * pi ^ 2)) * asoln2(θ, ω, Ω, D2)
    res = quadgk(f_int2, 0, 2 * pi, atol = errtol)[1]
    return res
end

ωn = 0.5
Ω0 = 3
D1 = 5
D2 = 8
testx = integral1(ωn, Ω0, D1, D2) + integral2(ωn, Ω0, D1, D2)

@time r1 = map(ωn -> integral1(ωn, Ω0, D1, D2) + integral2(ωn, Ω0, D1, D2), 0:0.01: 2 * pi)

plot(0:0.01: 2 * pi, real.(r1))

function Ξ_jl(z, system, j, l)
    Imp_j = system.Imps[j]
    Imp_l = system.Imps[l]
    m = system.M
    Ω0 = sqrt(2 * (system.k + system.α) / m)
    nmj = p2nm(system.N, Imp_j.p)
    nml = p2nm(system.N, Imp_l.p)
    D1 = nmj[1] - nml[1]
    D2 = nmj[2] - nml[2]

    return m * (integral1(- 1im * z, Ω0, D1, D2) + integral2(- 1im * z, Ω0, D1, D2))
end

@time Ξ_jl(1im * 0.5, testSys, 1, 1)

function Ξ(z, system)
    nImps = length(system.Imps)
    js = repeat(1:nImps, 1, nImps)
    ls = transpose(js)
    res = map((j, l) -> Ξ_jl(z, system, j, l), js, ls)
    return res
end

@time Ξ(1im * 0.5, testSys)

# Interaction energy (for now computed at T = 0)
# FI Integrand for zero-T energy
function F_I_Integrand_T0SL(z, system)
    nImps = length(system.Imps)
    Ξ_ = Ξ(z, system)
    λs = map(x -> x.λ, system.Imps)
    m = system.M
    ms = map(x -> m, system.Imps)
    α = Diagonal(1 ./ ms .- 1 ./ λs)
    Δ = (Matrix{Int}(I, nImps, nImps) .+ Ξ_ * α)
    Δ0_Inv = diag(Δ) |> Diagonal |> inv
    return (det(Δ0_Inv * Δ) |> Complex |> log |> real)
end

# 4e-4
@time F_I_Integrand_T0SL(1im * 0.5, testSys)

function F_I_T0(system)
    f(w) = F_I_Integrand_T0SL(1im * w, system)
    res = quadgk(f, 10, 0, 10, maxevals = NumEvals, rtol = ν)[1]
    return (res / (2 * π) / 2) # See the F_I formula for the division by 2
end

@time F_I_T0(testSys)

N = 20
k = 1
α = 2
m = 3
λ = 20
N_List = 2 : 3: 50
testImp1 = ImpuritySL(1, λ)
testImp2 = ImpuritySL(2, λ)
testSys = SystemSL(m, [testImp1, testImp2], k, α, T, N)

@time F_I_T0(testSys)


λs = map(x -> x.λ, system.Imps)
ms = map(x -> system.Ms[x.n], system.Imps)
α = Diagonal(1 ./ ms .- 1 ./ λs)


# Interaction energy (for now computed at T = 0)
# FI Integrand for zero-T energy
function F_I_Integrand_T0(z, system)
    nImps = length(system.Imps)
    Ξ_ = Ξ(z, system)
    λs = map(x -> x.λ, system.Imps)
    ms = map(x -> system.Ms[x.n], system.Imps)
    α = Diagonal(1 ./ ms .- 1 ./ λs)
    Δ = (Matrix{Int}(I, nImps, nImps) .+ Ξ_ * α)
    Δ0_Inv = diag(Δ) |> Diagonal |> inv
    return (det(Δ0_Inv * Δ) |> Complex |> log |> real)
end
