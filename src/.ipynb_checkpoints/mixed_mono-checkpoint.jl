using QuadGK
using Plots
using LaTeXStrings
using LinearAlgebra


# This file defines the function for the T = 0 computation for the
# mixed monoatomic chain

## Parameters
const ν = 1e-3;         # Relative tolerance for integration
const η = 1e-7;         # Small number used for i0
const NumEvals = 1e7;   # Maximum number of evaluations in quadgk
T = 1e-8

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


# Free Energy from Exact Diagonalization for two impurities
function Exact_Free_Energy(system)
    # easy refernce to system variables
    N = system.N
    m = system.m
    Ω_0 = system.Ω_0
    γ = system.γ
    λ = system.λ
    t = system.t
    T = system.T

    # Prepare a pristine chain potential energy matrix
    dv = 2 .* ones(N)
    ev = -ones(N - 1)

    U_Mat = SymTridiagonal(dv, ev) + zeros(N, N)
    # Make the system periodic by coupling the first and the last mass
    U_Mat[1,  N] = -1
    U_Mat[N,  1] = -1

    # add the harmonic perturbation
    U_Mat[t, t] = γ

    # Prepare the matrix of masses
    M_Mat = Diagonal(ones(N))
    # Replace ending pristine masses by impurities
    M_Mat[N, N] = 1 / λ

    # Calculate the eigenvalues of the system which give
    # the squares of the energies
    Ω2 = real(eigvals(inv(M_Mat) * U_Mat))
    # These Ω's are given in units of sqrt(K / m). We, however, want to give
    # the energies in the units of 2*sqrt(K / m). Hence, we divide the
    # energies by 2:
    Ω = real(sqrt.(complex(Ω2[2 : N])) / 2)
    # Ω = sqrt.(Ω2[2 : N]) / 2

    # The free energy for each mode consists of the vacuum portion Ω / 2 and
    # the finite-T portion T * log(1 - exp(-Ω / T))
    filter!(x -> x > 1e-12, Ω)
    total_energy = T * sum(log.(1 .- exp.(-Ω ./ T))) + sum(Ω) / 2

    return total_energy
end

function Exact_Int_Energy(system)

    # easy refernce to system variables
    N = system.N
    m = system.m
    Ω_0 = system.Ω_0
    γ = system.γ
    λ = system.λ
    t = system.t
    T = system.T

    F_I = Exact_Free_Energy(system)

    F_m = real(Exact_Free_Energy(Mixed_System_Single_Ptb(N, m, Ω_0, γ, 0, t, T)))
    F_h = real(Exact_Free_Energy(Mixed_System_Single_Ptb(N, m, Ω_0, 0, λ, t, T)))
    F_0 = real(Exact_Free_Energy(Mixed_System_Single_Ptb(N, m, Ω_0, 0, 0, t, T)))

    return (F_I - F_0) - (F_m - F_0) - (F_h - F_0)

end
