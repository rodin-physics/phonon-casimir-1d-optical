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

function H_0(Ms, N, Ω_0)
    m = Ms[1]
    M = Ms[2]
    H_0 = zeros(Float64, 2N, 2N)
    H_0[1, 1] = -2 / m
    H_0[1, 2] = 1 / m
    # H_0[1, 2N] = 1 / m # periodic boundary
    H_0[2N, 2N] = -2 / M
    H_0[2N, 2N - 1] = 1 / M
    # H_0[2N, 1] = 1 / M # periodic boundary
    for i in 2:(2N - 1)
        if i % 2 == 1
            H_0[i, i] = -2 / m
            H_0[i, i - 1] = 1 / m
            H_0[i, i + 1] = 1 / m
        elseif i % 2 == 0
            H_0[i, i] = -2 / M
            H_0[i, i - 1] = 1 / M
            H_0[i, i + 1] = 1 / M
        end
    end
    H_0 = Ω_0 * H_0
    return H_0
end

function H_I(Ms, N, Imps)
    nImps = length(Imps)
    H_I = zeros(Float64, 2N, 2N)
    # α = map(x -> 1 / Ms[x.n] - 1 / x.λ, Imps) #|> Diagonal
    for imp in Imps
        r = 2 * (imp.pos - 1) + imp.n
        αr = 1 / Ms[imp.n] - 1 / imp.λ
        if r == 1
            H_I[r, r] = 2 * αr
            H_I[r, r + 1] = - αr
            H_I[r, 2N] = - αr
        elseif r == 2N
            H_I[r, r] = 2 * αr
            H_I[r, r - 1] = - αr
            H_I[r, 1] = - αr
        else
            H_I[r, r] = 2 * αr
            H_I[r, r + 1] = - αr
            H_I[r, r - 1] = - αr
        end
    end
    return H_I
end

function testAndPlotRawDisp()
    N = 500
    testH0 = H_0([1, 2], N, 1)
    # # println(eigvals(testH0))
    # print(maximum(eigvals(testH0)))
    exactDispersion = map(x -> sqrt(abs(x)), eigvals(testH0))
    opticalbranch = reverse(sort(exactDispersion[1:N]))
    acousticbranch = sort(exactDispersion[N + 1: 2N])
    plot(1:N, acousticbranch)
    plot!(1:N, opticalbranch)
    #
    calcAcousticBranch = map(n ->  sqrt(Q2_sθ([1, 2], 0, pi * n / N)), 1:N)
    calcOptBranch = map(n ->  sqrt(Q2_sθ([1, 2], 1, pi * n / N)), 1:N)
    plot!(1:N, calcAcousticBranch)
    plot!(1:N, calcOptBranch)
end
##

testAndPlotRawDisp()
# periodic boundary conditions causes something weird to happen at the edges

function twoImpurities(i, j, s1, s2, impM)
    imps = Array{Any}(undef, 2)
    imps[1] = Impurity(i, s1, impM)
    imps[2] = Impurity(j, s2, impM)
    return imps
end



N = 500
testH0 = H_0([1, 1], N, 1)
testHI = H_I([1, 1], N,  twoImpurities(120, 130, 1, 1, 1))
println("\n")
println(testHI)

N = 500
testH0 = H_0([1, 1], N, 1)
testHI = H_I([1, 1], N,  twoImpurities(249, 251, 1, 1, 20))
testfullH = testH0 + testHI
exactDispersion = sort(map(x -> sqrt(-x), eigvals(testfullH)))
opticalbranch = exactDispersion[1:N]
acousticbranch = reverse(exactDispersion[N + 1: 2N])
plot(1:N, acousticbranch)
plot!(1:N, opticalbranch)

exactDispersion = sort(map(x -> sqrt(-x), eigvals(testH0)))
opticalbranch = exactDispersion[1:N]
acousticbranch = reverse(exactDispersion[N + 1: 2N])
plot!(1:N, acousticbranch)
plot!(1:N, opticalbranch)

### monoatomic ZPE test
function monoatZPE(N, m)
    H0 = H_0(Int(N / 2), 1, [m, m])
    exactZPEDispersion = map(x -> 0.5 * sqrt(-x), eigvals(H0))
    ZPE = sum(exactZPEDispersion)
    return ZPE
end

function monoatZPEw2Imp(N, m, i, j)
    H0 = H_0([m, m], Int(N / 2), 1)
    si = Int(i % 2)
    sj = Int(j % 2)
    truei = Int((i + si) / 2)
    truej = Int((j + sj) / 2)
    si = 2 - si
    sj = 2 - sj
    HI = H_I([m, m], Int(N / 2),  twoImpurities(truei, truej, si, sj, 1))
    # ZPE = 0.5*sum(map(x -> 0.5 * sqrt(abs(x)), eigvals(H0 + HI)))
    # ZPE = 0.5 * sum(eigvals(H0 + HI))
    # # ZPE = 0.5 * sum(eigvals(H0))
    # return ZPE
    # return sort(eigvals(H0 + HI))
    # return sort(map(x -> sqrt(-x), eigvals(H0 + HI)))
    ZPE = 0.5 * sum(map(x -> sqrt(-x), eigvals(H0 + HI)))
    return ZPE
end

function testAndPlot2impDisp()
    N = 500
    m = 1
    i = 249
    j = 251
    H0 = H_0([m, m], Int(N / 2), 1)
    si = Int(i % 2)
    sj = Int(j % 2)
    truei = Int((i + si) / 2)
    truej = Int((j + sj) / 2)
    si = 2 - si
    sj = 2 - sj
    HI = H_I([m, m], Int(N / 2),  twoImpurities(truei, truej, si, sj, 1))
    testH0 = H0 + HI
    # # println(eigvals(testH0))
    # print(maximum(eigvals(testH0)))
    exactDispersion = sort(map(x -> sqrt(abs(x)), eigvals(testH0)))
    opticalbranch = exactDispersion[1:Int(N / 2)]
    acousticbranch = reverse(exactDispersion[Int(N / 2) + 1: N])
    # plot(1:Int(N / 2), acousticbranch)
    # plot!(1:Int(N / 2), opticalbranch)
    # calcAcousticBranch = map(n ->  sqrt(Q2_sθ([1, 1], 0, pi * n / N)), 1:Int(N / 2))
    # calcOptBranch = map(n ->  sqrt(Q2_sθ([1, 1], 1, pi * n / N)), 1:Int(N / 2))
    # plot!(1:Int(N / 2), calcAcousticBranch)
    # plot!(1:Int(N / 2), calcOptBranch)
    return exactDispersion
end

testZPE = monoatZPEw2Imp(500, 1, 249, 251)
exactDispersion = monoatZPEw2Imp(500, 1, 249, 251)
opticalbranch = exactDispersion[1:250]
acousticbranch = reverse(exactDispersion[251: 500])
plot(1:250, acousticbranch)
plot!(1:250, opticalbranch)

exactDispersion = sort(map(x -> sqrt(-x), eigvals(H_0([1, 1], Int(500 / 2), 1))))
opticalbranch = exactDispersion[1:250]
acousticbranch = reverse(exactDispersion[251: 500])
plot!(1:250, acousticbranch)
plot!(1:250, opticalbranch)


function testMonoat2Imp()
    N = 1000
    i = 499
    j = 501
    distances = zeros(20)
    casimirEnergy = zeros(20)
    d = 2
    for n = 1:20
        distances[n] = d
        # println(typeof(monoatZPEw2Imp(N, 1, i, j)))
        casimirEnergy[n] = monoatZPEw2Imp(N, 1, i, j)
        i = i - 1
        j = j + 1
        d = d + 2
    end
    x = 1:20
    return casimirEnergy
end

y1 = testMonoat2Imp()
y1 = map(x -> x - maximum(y1), y1)
x = 1:20
plot(x, y1)
y1 = casimirAgainstDist() # base comp 400 600
y2 = casimirAgainstDist() # base comp 2 998
x = 1:20
plot(x, y)
py = map(v -> log(- v), y)
logx = map(v -> log(v), x)
# plot!(x, cb)
plot(logx, py)
