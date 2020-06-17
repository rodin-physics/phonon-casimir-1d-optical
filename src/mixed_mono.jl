using QuadGK
using Plots
using LaTeXStrings
using LinearAlgebra

## Parameters
const ν = 1e-3;         # Relative tolerance for integration
const η = 1e-7;         # Small number used for i0
const NumEvals = 1e7;   # Maximum number of evaluations in quadgk

function RtGnRS(N, m, z, Ω_0, γ, λ, t)
    # m is atom mass
    # z is matsubara freq
    # Ω_0 is dispersion constant
    # γ is mass perturbation
    # λ is potential well
    # t is location of potential well
    Asq = - ( z / Ω_0 ) ^ 2
    b = 2 * Asq + 1
    rminus = b - sqrt(b ^ 2 - 1)
    RtGnRS = zeros(Complex{Float64}, 2, 2)
    RtGnRS[1, 1] = - 2 * pi * m * γ * (2 - rminus - 1 / rminus)
    RtGnRS[1, 2] = - 8 * z * pi * λ / (Ω_0^2 * rminus ^ t)
    RtGnRS[2, 1] = - 8 * z * pi * λ * rminus ^ t / (Ω_0^2 )
    RtGnRS[2, 2] = - 4 * pi * λ / (m * Ω_0^2)
    RtGnRS = RtGnRS ./ (N * sqrt(b^2 - 1))
    return RtGnRS
end

function the_big_det(N, m, z, Ω_0, γ, λ, t)
    inside = Matrix{Complex{Float64}}(I, 2, 2) - RtGnRS(N, m, z, Ω_0, γ, λ, t)
    return det(inside)
end

function F_I_Integrand_T0(N, m, z, Ω_0, γ, λ, t, T)
    return log(the_big_det(N, m, z, Ω_0, γ, λ, t))
end

function F_I_T0(N, m, z, Ω_0, γ, λ, t, T)
    f(w) = F_I_Integrand_T0(N, m,  1im * w, Ω_0, γ, λ, t, T)
    # pole at 0

    res = quadgk(f, -Inf, 1e-7, Inf, maxevals = NumEvals, rtol = ν)[1]
    return (res / (2 * π) / 2) # See the F_I formula for the division by 2
end

N = 500
m = 1.
z = 1.5 * 1im
Ω_0 = 2.
γ = 3.
λ = 3.
t = 16

Asq = - ( z / Ω_0 ) ^ 2
b = 2 * Asq + 1
rminus = b - sqrt(b ^ 2 - 1)
- 8 * z * pi * λ / (Ω_0^2 * rminus ^ t)

@time RtGnRS(500, 1., 1.5 * 1im, 2., 3., 3., 2.)
print( RtGnRS(500, 1., 1.5 * 1im, 2., 3., 3., 2.))
typeof(z)
iden(Complex{Float64}, 2, 2)

@time the_big_det(N, m, z, Ω_0, γ, λ, t)
@time F_I_Integrand_T0(N, m, z, Ω_0, γ, λ, t, 1e-10)
@time F_I_T0(N, m, z, Ω_0, γ, λ, t, 1)

T = 1e-10
f(w) = F_I_Integrand_T0(N, m, 1im * w, Ω_0, γ, λ, t, T)

f(-1.0)
f(-0.99609375 + 0im)
res = quadgk(f, -Inf, 0, Inf, maxevals = NumEvals, rtol = ν)[1]

xVals = LinRange(-2.0, 2.0, 4000)
yVals = map(x -> real(f(x)), xVals)

pyplot()
plot(xVals, yVals)

f(0.0000001)
