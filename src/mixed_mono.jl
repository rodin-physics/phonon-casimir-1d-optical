using QuadGK
using Plots
using LaTeXStrings
using LinearAlgebra

## Parameters
const ν = 1e-3;         # Relative tolerance for integration
const η = 1e-7;         # Small number used for i0
const NumEvals = 1e7;   # Maximum number of evaluations in quadgk

struct Mixed_System_Single_Ptb
    N::Int
    m::Float64
    Ω_0::Float64
    γ::Float64
    λ::Float64
    t::Int # position of tweeze
    T::Float64
end

function F_I_T0_integrand(z, system)

    # easy refernce to system variables
    N = system.N
    m = system.m
    Ω_0 = system.Ω_0
    γ = system.γ
    λ = system.λ
    t = system.t

    # convenience variables
    Asq = - ( z / Ω_0 ) ^ 2
    b = 2 * Asq + 1
    sqbm1 = sqrt(b ^ 2 - 1)
    rmin = b - sqbm1

    # integrals done by hand
    I_plus = (2 * rmin ^ t) / (Ω_0^2 * sqbm1)
    I_minus = 2 / (rmin ^ t * Ω_0^2 * sqbm1)
    I_0 = 2 / (Ω_0^2 * sqbm1)
    I_1 = 1 - sqrt(Asq / (1 + Asq))
    I_2 = (rmin ^ t - rmin ^ (t + 1) / 2 - rmin ^ (t - 1) / 2) / sqbm1

    # the matrix
    F_mat = zeros(Complex{Float64}, 4, 4)
    # first column
    F_mat[1, 1] = 1 - γ * I_0
    F_mat[2, 1] = - γ * I_minus
    F_mat[3, 1] = - γ * z * I_0
    F_mat[4, 1] = - γ * z * I_minus
    # extraneous ones
    F_mat[2, 2] = 1
    F_mat[3, 3] = 1
    # last column
    F_mat[4, 1] = - λ * z * I_plus
    F_mat[4, 2] = - λ * z * I_0
    F_mat[4, 3] = - λ * I_2
    F_mat[4, 4] = 1 - λ * Ω_0^2 * I_1

    F_I_T0_integrand = log(det(F_mat))

    return F_I_T0_integrand
end

function F_I_T0(system)
    f(w) = F_I_T0_integrand(1im * (w + 1im * η), system)
    # pole at 0

    res = quadgk(f, -Inf, 0, Inf, maxevals = NumEvals, rtol = ν)[1]
    return (res / (2 * π) / 2) # See the F_I formula for the division by 2
end

function F_int_Integrand_T0(z, system)
    # easy refernce to system variables
    N = system.N
    m = system.m
    Ω_0 = system.Ω_0
    γ = system.γ
    λ = system.λ
    t = system.t

    # convenience variables
    Asq = - ( z / Ω_0 ) ^ 2
    b = 2 * Asq + 1
    sqbm1 = sqrt(b ^ 2 - 1)
    rmin = b - sqbm1

    # integrals done by hand
    I_plus = (2 * rmin ^ t) / (Ω_0^2 * sqbm1)
    I_minus = 2 / (rmin ^ t * Ω_0^2 * sqbm1)
    I_0 = 2 / (Ω_0^2 * sqbm1)

    sub_mat_1a = zeros(Complex{Float64}, 2, 2)
    sub_mat_1a[1, 1] = 1 / m
    sub_mat_1a[1, 2] = z
    sub_mat_1a[2, 1] = z
    sub_mat_1a[2, 2] = m * z ^ 2

    sub_mat_1b = zeros(Complex{Float64}, 4, 4)
    sub_mat_1b[3, 3] = m
    sub_mat_1b[4, 4] = m

    Y_D = zeros(Complex{Float64}, 2, 2)
    Y_D[1, 1] = I_0
    Y_D[2, 2] = I_0

    Y_off = zeros(Complex{Float64}, 2, 2)
    Y_off[1, 2] = I_plus
    Y_off[2, 1] = I_minus

    int_mat = zeros(Complex{Float64}, 4, 4)
    int_mat[1, 1] = γ
    int_mat[4, 4] = λ

    Ξ_0n = kron(sub_mat_1a, Y_D) + sub_mat_1b
    Ξ_In = kron(sub_mat_1a, Y_off)

    inner_arg = I + inv(I + Ξ_0n * int_mat) * Ξ_In * int_mat
    integrand = log(det(inner_arg))

    return integrand

end

function F_int_T0(system)
    f(w) = F_int_Integrand_T0(1im * (w + 1im * η), system)
    # pole at 0

    res = quadgk(f, -Inf, 0, Inf, maxevals = NumEvals, rtol = ν)[1]
    return (res / (2 * π) / 2) # See the F_I formula for the division by 2
end

N = 500
m = 1.
Ω_0 = 4.
γ = 0.3
λ = 0.3
T = 1e-10
t = 16

pyplot()

test_sys = Mixed_System_Single_Ptb(N, m, Ω_0, γ, λ, t, T)
# 1e-5s
z = 1.5 * 1im
@time F_I_T0_integrand(z, test_sys)
# 1e-3s
@time F_I_T0(test_sys)

# 1e-4s
@time F_int_Integrand_T0(z, test_sys)
# 6e-4s
@time F_int_T0(test_sys)
# 0.09s
@time range_of_F = map(t -> real(F_int_T0(Mixed_System_Single_Ptb(N, m, Ω_0, γ, λ, t, T))), 2:50)
# evidently some spikes
plot(2:50, range_of_F)



# 0.09s
@time range_of_F = map(t -> real(F_int_T0(Mixed_System_Single_Ptb(N, m, Ω_0, γ, λ, t, T))), 50:74)
# evidently some spikes
plot(50:74, range_of_F)

# 0.09s
@time range_of_F = map(t -> real(F_int_T0(Mixed_System_Single_Ptb(N, m, Ω_0, γ, λ, t, T))), 100:120)
# evidently some spikes
plot(50:74, range_of_F)


x = 4
# 1e-1s
@time range_of_F = map(t -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, Ω_0, γ, λ, t, T))), 2:50)
# evidently no dependence on t
plot(2:50, range_of_F)

# check if any interaction energy between mass and harmonic perturbations
Fint = real(F_I_T0(Mixed_System_Single_Ptb(N, m, Ω_0, γ, λ, t))
              - F_I_T0(Mixed_System_Single_Ptb(N, m, Ω_0, γ, 0, t))
              - F_I_T0(Mixed_System_Single_Ptb(N, m, Ω_0, 0, λ, t)))

λs = [0.00001, 0.0001, 0.001, 0.01, 0.1, 10, 100.]
γs = [0.00001, 0.0001, 0.001, 0.01, 0.1, 10, 100.]

# γ dependence for Ω_0 = 0.01
pyplot()
range_of_F_1 = map(λ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 0.01, γs[1], λ, t))), λs)
range_of_F_3 = map(λ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 0.01, γs[3], λ, t))), λs)
range_of_F_5 = map(λ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 0.01, γs[5], λ, t))), λs)
range_of_F_7 = map(λ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 0.01, γs[7], λ, t))), λs)
plot(λs, range_of_F_1)
plot!(λs, range_of_F_3)
plot!(λs, range_of_F_5)
plot!(λs, range_of_F_7)

# γ dependence for Ω_0 = 1
pyplot()
range_of_F_1 = map(λ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 1., γs[1], λ, t))), λs)
range_of_F_3 = map(λ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 1., γs[3], λ, t))), λs)
range_of_F_5 = map(λ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 1., γs[5], λ, t))), λs)
range_of_F_7 = map(λ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 1., γs[7], λ, t))), λs)
plot(λs, range_of_F_1)
plot!(λs, range_of_F_3)
plot!(λs, range_of_F_5)
plot!(λs, range_of_F_7)

# γ dependence for Ω_0 = 100
pyplot()
range_of_F_1 = map(λ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 100., γs[1], λ, t))), λs)
range_of_F_3 = map(λ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 100., γs[3], λ, t))), λs)
range_of_F_5 = map(λ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 100., γs[5], λ, t))), λs)
range_of_F_7 = map(λ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 100., γs[7], λ, t))), λs)
plot(λs, range_of_F_1)
plot!(λs, range_of_F_3)
plot!(λs, range_of_F_5)
plot!(λs, range_of_F_7)

# λ dependence for Ω_0 = 0.01
pyplot()
range_of_F_1 = map(γ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 0.01, γ, λs[1], t))), γs)
range_of_F_3 = map(γ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 0.01, γ, λs[3], t))), γs)
range_of_F_5 = map(γ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 0.01, γ, λs[5], t))), γs)
range_of_F_7 = map(γ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 0.01, γ, λs[7], t))), γs)
plot(γs, range_of_F_1)
plot!(γs, range_of_F_3)
plot!(γs, range_of_F_5)
plot!(γs, range_of_F_7)

# λ dependence for Ω_0 = 1.
pyplot()
range_of_F_1 = map(γ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 1., γ, λs[1], t))), γs)
range_of_F_3 = map(γ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 1., γ, λs[3], t))), γs)
range_of_F_5 = map(γ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 1., γ, λs[5], t))), γs)
range_of_F_7 = map(γ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 1., γ, λs[7], t))), γs)
plot(γs, range_of_F_1)
plot!(γs, range_of_F_3)
plot!(γs, range_of_F_5)
plot!(γs, range_of_F_7)

# λ dependence for Ω_0 = 100.
pyplot()
range_of_F_1 = map(γ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 100., γ, λs[1], t))), γs)
range_of_F_3 = map(γ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 100., γ, λs[3], t))), γs)
range_of_F_5 = map(γ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 100., γ, λs[5], t))), γs)
range_of_F_7 = map(γ -> real(F_I_T0(Mixed_System_Single_Ptb(N, m, 100., γ, λs[7], t))), γs)
plot(γs, range_of_F_1)
plot!(γs, range_of_F_3)
plot!(γs, range_of_F_5)
plot!(γs, range_of_F_7)

test = 4

# function RtGnRS(N, m, z, Ω_0, γ, λ, t)
#     # m is atom mass
#     # z is matsubara freq
#     # Ω_0 is dispersion constant
#     # γ is mass perturbation
#     # λ is potential well
#     # t is location of potential well
#     Asq = - ( z / Ω_0 ) ^ 2
#     b = 2 * Asq + 1
#     rminus = b - sqrt(b ^ 2 - 1)
#     RtGnRS = zeros(Complex{Float64}, 2, 2)
#     RtGnRS[1, 1] = - m * γ * (1 - sqrt(Asq / (1 + Asq)))
#     RtGnRS[1, 2] = - 2 * z * λ / (Ω_0^2 * rminus ^ t * sqrt(b ^ 2 - 1))
#     RtGnRS[2, 1] = - 2 * z * λ * rminus ^ t / (Ω_0^2 * sqrt(b ^ 2 - 1))
#     RtGnRS[2, 2] = - 2 * λ / (m * Ω_0^2 * sqrt(b ^ 2 - 1))
#     RtGnRS = RtGnRS ./ (N * sqrt(b^2 - 1))
#     return RtGnRS
# end
#
# function the_big_det(N, m, z, Ω_0, γ, λ, t)
#     inside = Matrix{Complex{Float64}}(I, 2, 2) - RtGnRS(N, m, z, Ω_0, γ, λ, t)
#     # Asq = - ( z / Ω_0 ) ^ 2
#     # inside = 1 + m * γ * (1 - sqrt(Asq / (1 + Asq)))
#     return det(inside)
# end
#
# function F_I_Integrand_T0(N, m, z, Ω_0, γ, λ, t, T)
#     return log(the_big_det(N, m, z, Ω_0, γ, λ, t))
# end
#
# function F_I_T0(N, m, z, Ω_0, γ, λ, t, T)
#     f(w) = F_I_Integrand_T0(N, m,  1im * (w + 1im * η), Ω_0, γ, λ, t, T)
#     # pole at 0
#
#     res = quadgk(f, -Inf, 0, Inf, maxevals = NumEvals, rtol = ν)[1]
#     return (res / (2 * π) / 2) # See the F_I formula for the division by 2
# end
#
# N = 500
# m = 1.
# z = 1.5 * 1im
# Ω_0 = 4.
# γ = 0.3
# λ = 0.3
# T = 1e-10
# t = 16
#
# @time the_big_det(N, m, z, Ω_0, γ, λ, t)
# @time F_I_Integrand_T0(N, m, 1.2 * 1im, Ω_0, γ, λ, 5, 1e-10)
# @time F_I_T0(N, m, z, Ω_0, γ, λ, 100, 1)
#
# f(t) = F_I_T0(N, m, z, Ω_0, γ, λ, t, 1)
#
# range_of_F = map(t -> real(f(t)), 1:3:31)
# plot(1:3:31, range_of_F)
# print(range_of_F)
#
# det_X = LinRange(-5e4, 5e4, 4000)
#
# det_Y = map(w -> real(the_big_det(N, m, 1im * w, Ω_0, γ, λ, t)), det_X)
# plot(det_X, det_Y)
#
#
# f(w) = F_I_Integrand_T0(N, m, 1im * (w + 1im * η), Ω_0, 0, λ, t, T)
#
# f(-1.0)
# f(-0.99609375 + 0im)
# @time res = quadgk(f, -1e5, 0, 1e5, maxevals = NumEvals, rtol = ν * 1e2, atol = 1e-9)[1]
#
# xVals = LinRange(-5e4, 5e4, 4000)
# yVals = map(x -> imag(f(x)), xVals)
#
# xIntBounds = LinRange(1e1, 1e4, 100)
# yIntValsReal = map(x-> real(quadgk(f, -x, 0, x, maxevals = NumEvals, rtol = ν, atol = 1e-9)[1]), xIntBounds)
# # yIntValsReal = map(x -> real(x), yIntVals)
# pyplot()
# plot(xVals, yVals)
# plot(xIntBounds, yIntValsReal)
#
# f(0.000000001 + 1im * η)
# f(0.000000001)
# f(4000)
