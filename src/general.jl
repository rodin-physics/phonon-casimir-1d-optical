using Plots
using QuadGK
using LinearAlgebra

## Parameters
const ν = 1e-3;         # Relative tolerance for integration
const η = 1e-7;         # Small number used for i0
const NumEvals = 1e7;   # Maximum number of evaluations in quadgk

## Functions
function Q2_sθ(m, M, s, θ)
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
function e_sθ(m, M, s, θ)
    return ([1 + exp(-1im * θ), 2 - m * Q2_sθ(m, M, s, θ)])
end

# Mode couling term
function coupling(m, M, s, θ)
    M_mat = [m 0; 0 M]
    e = e_sθ(m, M, s, θ)
    return (M_mat * e * e' * M_mat) / (e' * M_mat * e)
end

# Impurity coupling
function YGY(m, M, z, D, Aj, Ak)
    f_int_acc(θ) =
        Q2_sθ(m, M, 0, θ) * exp(1im * D * θ) / (-z^2 + Q2_sθ(m, M, 0, θ)) *
        (Aj' * coupling(m, M, 0, θ) * Ak)

    f_int_opt(θ) =
        Q2_sθ(m, M, 1, θ) * exp(1im * D * θ) / (-z^2 + Q2_sθ(m, M, 1, θ)) *
        (Aj' * coupling(m, M, 1, θ) * Ak)

    f_int(θ) = f_int_opt(θ) + f_int_acc(θ)
    res = quadgk(f_int, 0, 2 * π)[1]
    return (-res) / (2 * π)
end
