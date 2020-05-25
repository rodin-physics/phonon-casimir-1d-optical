using LinearAlgebra
using Plots
pyplot()

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

get_adj(3, 1, 'c', -1)



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

# just a test function to check that exact diagonalisation works without impurities
function SLED(N, k, α, m)
    Vx = zeros(N ^ 2, N ^ 2)
    Vy = zeros(N ^ 2, N ^ 2)
    M_Mat = fill(m, N ^ 2) |> Diagonal
    for p in 1: N ^ 2
        # println(p)
        Vx[p, p] = 2 * k + 2 * α
        left = get_adj(N, p, 'r', -1)
        right = get_adj(N, p, 'r', 1)
        Vx[p, left] = - k
        Vx[p, right] = - k
        top = get_adj(N, p, 'c', -1)
        bottom = get_adj(N, p, 'c', 1)
        Vx[p, top] = - α
        Vx[p, bottom] = -α
    end
    for p in 1: N ^ 2
        # println(p)
        Vy[p, p] = 2 * k + 2 * α
        left = get_adj(N, p, 'r', -1)
        right = get_adj(N, p, 'r', 1)
        Vy[p, left] = - α
        Vy[p, right] = - α
        top = get_adj(N, p, 'c', -1)
        bottom = get_adj(N, p, 'c', 1)
        Vy[p, top] = - k
        Vy[p, bottom] = - k
    end
    ω_x = sqrt.(abs.(eigvals(inv(M_Mat) * Vx)))
    ω_y = sqrt.(abs.(eigvals(inv(M_Mat) * Vy)))
    return(0.5 * sum(ω_x) + 0.5 * sum(ω_y))
end

function SLED_Imp(system)
    M = system.M
    Imps = system.Imps
    k = system.k
    α = system.α
    T = system.T
    N = system.N

    Vx = zeros(N ^ 2, N ^ 2)
    Vy = zeros(N ^ 2, N ^ 2)
    M_Mat = fill(M, N ^ 2)
    for ii in Imps
        coord = ii.p
        M_Mat[coord] = ii.λ
    end
    M_Mat = Diagonal(M_Mat)
    for p in 1: N ^ 2
        # println(p)
        Vx[p, p] = 2 * k + 2 * α
        left = get_adj(N, p, 'r', -1)
        right = get_adj(N, p, 'r', 1)
        Vx[p, left] = - k
        Vx[p, right] = - k
        top = get_adj(N, p, 'c', -1)
        bottom = get_adj(N, p, 'c', 1)
        Vx[p, top] = - α
        Vx[p, bottom] = -α
    end
    for p in 1: N ^ 2
        # println(p)
        Vy[p, p] = 2 * k + 2 * α
        left = get_adj(N, p, 'r', -1)
        right = get_adj(N, p, 'r', 1)
        Vy[p, left] = - α
        Vy[p, right] = - α
        top = get_adj(N, p, 'c', -1)
        bottom = get_adj(N, p, 'c', 1)
        Vy[p, top] = - k
        Vy[p, bottom] = - k
    end
    ω_x = sqrt.(abs.(eigvals(inv(M_Mat) * Vx)))
    ω_y = sqrt.(abs.(eigvals(inv(M_Mat) * Vy)))
    filter!(x -> x > 1e-12, ω_x)
    filter!(x -> x > 1e-12, ω_y)
    U = 0.5 * sum(ω_x) + 0.5 * sum(ω_y)
    total_energy = U + T * (sum(log.(1 .- exp.(-ω_x ./ T))) + sum(log.(1 .- exp.(-ω_y ./ T))))
    return total_energy
end

function SLZPE_0(system)
    M = system.M
    k = system.k
    α = system.α
    T = system.T
    N = system.N

    dispersion_x = zeros(N ^ 2)
    dispersion_y = zeros(N ^ 2)
    for qx in 1:N
        for qy in 1:N
            ω_x = sqrt((2 * k / m) * (1 - cos(qx)) + (2 * α / m) * (1 - cos(qy)))
            dispersion_x[nm2p(N, qx, qy)] = ω_x
            ω_y = sqrt((2 * k / m) * (1 - cos(qy)) + (2 * α / m) * (1 - cos(qx)))
            dispersion_y[nm2p(N, qx, qy)] = ω_y
        end
    end

    U = 0.5 * sum(dispersion_x) + 0.5 * sum(dispersion_y)
    total_energy = U + T * (sum(log.(1 .- exp.(-dispersion_x ./ T))) + sum(log.(1 .- exp.(-dispersion_y ./ T))))
    return total_energy
end

function casimirEDSL(system)
    U_0 = SLZPE_0(system)
    U_imp = SLED_Imp(system)

    U_casimir = U_imp - U_0
    for imp in system.Imps
        system_i = SystemSL(system.M, [imp], system.k, system.α, system.T, system.N)
        U_i = SLED_Imp(system_i) - U_0
        U_casimir -= U_i
    end

    return U_casimir
end

N = 30
k = 1
α = 2
m = 3
λ = 20
N_List = 2 : 3: 50
T = 1e-12;
testImp1 = ImpuritySL(1, λ)
testImp2 = ImpuritySL(5, λ)
testSys = SystemSL(m, [testImp1, testImp2], k, α, T, N)


SLZPE_0(testSys)
SLED_Imp(testSys)

@time begin
casimirEDSL(testSys)
end

row_crawl = 2: 1: 10
diag_crawl = N + 2: N + 2: 10 * (N + 2)

@time begin
r1 = map(d -> casimirEDSL(SystemSL(m, [ImpuritySL(1, λ), ImpuritySL(d, λ)], k, α, T, N)), row_crawl)
end

@time begin
r2 = map(d -> casimirEDSL(SystemSL(m, [ImpuritySL(1, λ), ImpuritySL(d, λ)], k, α, T, N)), diag_crawl)
end

diagDist = map(d -> d * sqrt(2), 1:10)
plot(row_crawl, r1)
plot!(diagDist, r2)

plot(N_List, r1)
plot!(N_List, r2)

testHx = squareLatticeH(N, k, α)

eigsHx = eigvals(testHx)
eigsHx[1, 1] = 0
plot(1:N^2, sqrt.(eigsHx))

analyticDisp_x, analyticDisp_y = SLanalyticDisp(N, k, α)
plot(1:N^2, analyticDisp_x)
plot!(1:N^2, analyticDisp_y)




plot([1, 2, 3, 4, 5], [1, 2, 3, 4, 5])

pyplot()
plot([1, 2, 3, 4, 5], [1, 2, 3, 4, 5])
