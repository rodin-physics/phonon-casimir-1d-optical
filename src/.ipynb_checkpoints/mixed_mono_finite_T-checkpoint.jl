using QuadGK
using Plots
using LaTeXStrings
using LinearAlgebra

## Parameters
const ν = 1e-1;         # Small number for relative tolerance
const α = 1e-3;         # Small number for absolute tolerance
const η = 1e-5;         # Small number for moving the contour off the real axis
const NumEvals = 1e7;   # Maximum number of evaluations in quadgk
# Impurity type
struct Impurity
    pos::Int        # Position of the impurity
    M::Float64      # Impurity mass in units of m
    Δ::Float64      # Harmonic potential perturbation in units of 4K
end
function propagator(z, D)
    x = z^2
    res = ((1 - 2 * x - 2 * √(-x) * √(1 - x))^D) / (√(-x) * √(1 - x))
    return res
end
function propagatorprime(z, D)
    if D > 0
        x1 = 1 - 2 * z^2 - 2 * sqrt( - z^2) * sqrt(1 - z^2)
        x2 = sqrt( - z^2) * ( 1 - 2 * z^2) * sqrt(1 - z^2) + 2 * D * z^2 * (-1 + z^2)
        x3 = z^3 * (-1 + z^2)^2
        upper = x1 ^ D * x2
        lower = x3
    elseif D == 0
        upper = z * (1 - 2 * z ^ 2)
        lower = (-z^2) ^ (3 / 2) * (1 - z^2 ) ^ (3/2)
    end
    return upper / lower
end

function Az(z, imps)
    nImps = length(imps)
    Λ_ = 1 ./ map(x -> x.M, imps) .- 1 |> Diagonal
    Δ_ = map(x -> x.Δ, imps) |> Diagonal
    Imp_Mat = repeat(imps, 1, nImps)    # Impurity matrix
    ImpT_Mat = permutedims(Imp_Mat)     # Transpose of impurity matrix
    Y = map((x, y) -> propagator(z, abs.(x.pos - y.pos)), Imp_Mat, ImpT_Mat)
    return [
        Y * Δ_ z .* Y * Λ_
        z .* Y * Δ_ (Matrix{Int}(I, nImps, nImps) + z^2 .* Y) * Λ_
    ]
end
function Azprime(z, imps)
    nImps = length(imps)
    Λ_ = 1 ./ map(x -> x.M, imps) .- 1 |> Diagonal
    Δ_ = map(x -> x.Δ, imps) |> Diagonal
    Imp_Mat = repeat(imps, 1, nImps)    # Impurity matrix
    ImpT_Mat = permutedims(Imp_Mat)     # Transpose of impurity matrix
    Y = map((x, y) -> propagator(z, abs.(x.pos - y.pos)), Imp_Mat, ImpT_Mat)
    dY = map((x, y) -> propagatorprime(z, abs.(x.pos - y.pos)), Imp_Mat, ImpT_Mat)
    return [
        dY * Δ_ (Y + z .* dY) * Λ_
        (Y + z .* dY) * Δ_ (2 * z .* Y + z^2 .* dY) * Λ_
    ]
end

function F_I_Integrand_T(z, imps, T)
    # print(z)
    nImps = length(imps)
    Az_mat = Az(z, imps)
    dAz_mat = Azprime(z, imps)
    term1 = log(exp(- z / T))
    term2 = tr(inv(Matrix{Int}(I, 2 * nImps, 2 * nImps) + Az_mat) * dAz_mat)
    return term1 * term2
end

function F_I(imps, T)
    f_plus(x) = F_I_Integrand_T(x + η * 1im, imps, T)
    f_min(x) = F_I_Integrand_T(x - η * 1im, imps, T)
    f(x) = f_plus(x) + f_min(x)
    res = quadgk(f, -Inf, Inf, maxevals = NumEvals)[1]
    return (res / (2 * π))
end

function F_I_int(t, imp_set, T)
    imp_m = Impurity(0, imp_set[1], 0)
    imp_t = Impurity(t, 1., imp_set[2])
    # imp_t1 = Impurity(1, 1., imp_set[2])
    # F1 = real(F_I([imp_m, imp_t1]) - F_I([imp_m]) - F_I([imp_t1]))
    F1 = 1.
    # F_I([imp_m, imp_t], T)
    # F_I([imp_m], T)
    # F_I([imp_t], T)
    return real(F_I([imp_m, imp_t], T) - F_I([imp_m], T) - F_I([imp_t], T)) / F1
end

function Exact_Free_Energy(N, imps, T)
    # Prepare a pristine chain potential energy matrix
    dv = 2 .* ones(N)
    ev = -ones(N - 1)
    U_Mat = SymTridiagonal(dv, ev) + zeros(N, N)
    # Make the system periodic by coupling the first and the last mass
    U_Mat[1, N] = -1
    U_Mat[N, 1] = -1
    # Add the external harmonic potential. The factor of 4 is because Δ is
    # given in units of 4 * K
    for ii = 1:length(imps)
        imp = imps[ii]
        U_Mat[imp.pos+1, imp.pos+1] += 4 * imp.Δ
    end
    # Prepare the matrix of masses
    M_Mat = Diagonal(ones(N))
    # Modify the masses
    for ii = 1:length(imps)
        imp = imps[ii]
        M_Mat[imp.pos+1, imp.pos+1] = imp.M
    end
    # Calculate the eigenvalues of the system which give
    # the squares of the energies
    Ω2 = real(eigvals(inv(M_Mat) * U_Mat))
    # These Ω's are given in units of sqrt(K / m). We, however, want to give
    # the energies in the units of 2*sqrt(K / m). Hence, we divide the
    # energies by 2:
    Ω = real(sqrt.(Ω2[2:N]) / 2)

    filter!(x -> x > 1e-12, Ω)
    total_energy = T * sum(log.(1 .- exp.(-Ω ./  T))) + sum(Ω) / 2

    return total_energy
end

function Exact_Int_Energy(N, t, imp_set, T)
    imp_m = Impurity(0, imp_set[1], 0)
    imp_t = Impurity(t, 1., imp_set[2])
    F_0 = Exact_Free_Energy(N, [], T)
    F_2 = Exact_Free_Energy(N, [imp_m, imp_t], T) - F_0
    F_m = Exact_Free_Energy(N, [imp_m], T) - F_0
    F_t = Exact_Free_Energy(N, [imp_t], T) - F_0
    return F_2 - F_m - F_t
end

function Exact_2m_Halfway(N, t, imp_set, T)
    imp_m1 = Impurity(0, imp_set[1], 0)
    imp_m2 = Impurity(t, imp_set[2], 0)
    imp_mH = Impurity(Int(N / 2), imp_set[2], 0)
    F_Halfway = Exact_Free_Energy(N, [imp_m1, imp_mH], T)
    F_Abs = Exact_Free_Energy(N, [imp_m1, imp_m2], T)
    return F_Abs - F_Halfway
end


t = 5
T = 0.1
imp_m = Impurity(0, 1.2, 0)
imp_t = Impurity(t, 1., 1.44)
@time F_I([imp_m, imp_t], T)
F_I_Integrand_T(0.9im + η, [imp_m, imp_t], T)
F_I_Integrand_T(-0.5 + η * 1im, [imp_m, imp_t], T)
F_I_Integrand_T(-0.5 - η * 1im, [imp_m, imp_t], T)
plot(-2:0.02:2, map(x -> real.(F_I_Integrand_T(x + η * 1im, [imp_m, imp_t], T)), -2:0.02:2))
# plot!(-2:0.1:2, map(x -> imag.(F_I_Integrand_T(x + η * 1im, [imp_m, imp_t], T)), -2:0.1:2))
Ds = 1:50
# T = 1.2
set_1a = [1.2, 1.44]
impset_1a = [imp_m, imp_t]
@time F_I([imp_m, imp_t], T)
@time F_I([imp_m], T)
@time F_I([imp_t], T)
@time F_I([imp_m, imp_t], T * 100)
real(F_I([imp_m, imp_t], T) - F_I([imp_m], T) - F_I([imp_t], T))
@time F_I_int(t, set_1a, T)
N = 500
@time Exact_Int_Energy(N, t, set_1a, T)

Ds_Exact = 2 : 3 : 20
F_I_int(10, set_1a, 1e-12)
T = 4
# range_1aM_exact = map(t -> Exact_2m_Halfway(N, t, set_1a, T), Ds_Exact)
# imp_m1 = Impurity(0, 1.2, 0)
# imp_m2 = Impurity(0, 1.44, 0)
range_1aM = map(t -> real(F_I([imp_m1, imp_m2], T) - F_I([imp_m1], T) - F_I([imp_m2], T)), Ds_Exact)

range_1a = map(t -> F_I_int(t, set_1a, T), Ds_Exact)
range_1a_exact = map(t -> Exact_Int_Energy(N, t, set_1a, T), Ds_Exact)

# plot(Ds_Exact, range_1aM)
# plot!(Ds_Exact, range_1aM_exact)

plot(Ds_Exact, range_1a)
plot!(Ds_Exact, range_1a_exact)

range_1a_exact_0T = map(t -> Exact_Int_Energy(N, t, set_1a, 1e-1), Ds_Exact)
plot!(Ds_Exact, range_1a_exact_0T)
