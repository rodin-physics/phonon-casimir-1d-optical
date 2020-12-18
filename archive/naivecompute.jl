using Plots
using QuadGK
using LinearAlgebra

## Parameters
const ν = 1e-3;         # Relative tolerance for integration
const η = 1e-7;         # Small number used for i0
const NumEvals = 1e7;   # Maximum number of evaluations in quadgk


Ω = 2
m = 1
M = 1
impW = 20
N = 500

# dispersion relation
function Ω_sq(m, M, s, q, Ω)
    # s should be 1 for optical branch, 0 for acoustic branch
    # - eventually generalise to n branches
    sqTerm = sqrt(m^2 + M^2 + 2 * m * M * cos(q)) / (m * M)
    invariantTerm = (m + M) / (m * M)
    if s == 0
        return Ω * sqrt(invariantTerm - sqTerm)
    elseif s == 1
        return Ω * sqrt(invariantTerm + sqTerm)
    else
        println("This branch is not available at the moment")
    end
end

# create an unnormalised basis vector
function basis_sq(m, M, s, q, Ω)
    esq1 = 1 + exp(1im * q )
    esq2 = 2 - Ω_sq(m, M, s, q, Ω) ^ 2
    return [esq1; esq2]
end

function YGY(m, M, j, k, z, N, numPts)
    # even j, k we set to be small masses,
    # odd j, k we set to be large masses
    d = j - k
    if j % 2 == 0
        Aj = [1; 0]
    elseif j % 2 == 1
        Aj = [0; 1]
    else
        println("Something wrong with j")
    end
    if k % 2 == 0
        Ak = [1; 0]
    elseif k % 2 == 1
        Ak = [0; 1]
    else
        println("Something wrong with j")
    end
    total = 0
    qaxis = range(-pi,stop= pi,length=N)
    mMat = [m 0; 0 M]
    s = 0
    for i in 1:N
        q = qaxis[i]
        Osq = Ω_sq(m,M,s,q,Ω)
        esq = basis_sq(m, M, 0, q, Ω)
        esqdag = conj(esq')
        factor0 = - 1 / N
        factor1 = exp(1im * d * q)
        factor2 = 1 / (esqdag * mMat * esq)[1, 1]
        factor3 = Aj' * mMat * esq * esqdag * mMat * Ak
        factor4 = Osq ^ 2 / (Osq ^ 2 - z^2) # IS IT MINUS Z ^ 2
        total = total + factor0 * factor1 * factor2 * factor3 * factor4
    end
    s = 1
    for i in 1:N
        q = qaxis[i]
        Osq = Ω_sq(m,M,s,q,Ω)
        esq = basis_sq(m, M, 0, q, Ω)
        esqdag = conj(esq')
        factor0 = - 1 / N
        factor1 = exp(1im * d * q)
        factor2 = 1 / (esqdag * mMat * esq)[1, 1]
        factor3 = Aj' * mMat * esq * esqdag * mMat * Ak
        factor4 = Osq ^ 2 / (Osq ^ 2 - z^2) # IS IT MINUS Z ^ 2
        total = total + factor0 * factor1 * factor2 * factor3 * factor4
    end
    return total
end

function twoFimpIntegrand(z, Ω, m, M, impM, j, k, N)
    alpha = [1/impM - 1/m 0; 0 1/impM - 1/m]
    YGYfactor = YGY(m, M, j, k, z, N, 10)
    integrand = log(det([1 0; 0 1] - YGYfactor * alpha))
    return integrand
end

function f(z)
    return twoFimpIntegrand(z, Ω, m, M, 0.7, 4, 6, 10)
end

function computeManyDistances()
    numD = 50
    j = 198
    k = 202
    fVals = zeros(Complex, numD)
    for i = 1:numD
        f(z) = twoFimpIntegrand(z, Ω, m, M, impW, j, k, N)
        FI = (1 / (4 * pi * 1im)) * quadgk(f, 0, 100000 / 1im)[1]
        # println("hi!")
        # println(FI)
        fVals[i] = FI
        j = j - 2
        k = k + 2
    end
    return fVals
end

Fimps = computeManyDistances()
yaxis = real(Fimps)

xaxis = range(0,stop= 50,length=50)
plot(xaxis, yaxis)

println(real(computeManyDistances()))

FI = (1 / (4 * pi * 1im)) * quadgk(f, 0, 10)[1]
println(FI)

println(testFimp())
println(f(5))

println(quadgk(f, 0, 0.5))


zaxis = range(0,stop= 100,length=numPts)
yplot = zeros(numPts)
for i in 1:numPts
    z = zaxis[i]
    yplot[i] = f(z)
end
plot(zaxis,yplot)
