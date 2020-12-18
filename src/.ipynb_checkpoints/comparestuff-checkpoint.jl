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
        (x, y) -> YGY(Ms, z, abs.(x.pos - y.pos), x.n, y.n),
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
YGY([1, 2], 1.2 + 1im * η, 2, 1,2)
# Q2_sθ([1,2.], 1, 1.2)

function H_0(N, Ω_0, Ms)
    m = Ms[1]
    M = Ms[2]
    H_0 = zeros(Float64, 2N, 2N)
    H_0[1, 1] = -2 / m
    H_0[1, 2] = 1 / m
    H_0[N, N] = -2 / M
    H_0[N, N - 1] = 1 / M
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

function twoImpurities(i, j, s1, s2, impM)
    imps = Array{Any}(undef, 2)
    imps[1] = Impurity(i, s1, impM)
    imps[2] = Impurity(j, s2, impM)
    return imps
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
        elseif r == 2 * N
            H_I[r, r] = 2 * αr
            H_I[r, r - 1] = - αr
        else
            H_I[r, r] = 2 * αr
            H_I[r, r + 1] = - αr
            H_I[r, r - 1] = - αr
        end
    end
    return H_I
end

function monoatZPEwImp(N, m, i, j)
    H0 = H_0(Int(N / 2), 1, [m, m])
    si = Int(i % 2)
    sj = Int(j % 2)
    truei = Int((i + si) / 2)
    truej = Int((j + sj) / 2)
    if si == 0
        si = 2
    end
    if sj == 0
        sj = 2
    end
    # println(truei)
    # println(truej)
    # println(si)
    # println(sj)
    HI = H_I([m, m], Int(N / 2),  twoImpurities(truei, truej, si, sj, 20))
    exactZPEDispersion = map(x -> 0.5 * sqrt(-x), eigvals(H0 + HI))
    ZPE = sum(exactZPEDispersion)
    return ZPE
end

monoatZPEwImp(1000, 1, 499, 501)
YGY([1, 1], 1.2 + 1im * η, 2, 1,2)
FI_Integrand([1, 1], 1.2 + 1im * η, twoImpurities(250, 251, 1, 1, 20))
