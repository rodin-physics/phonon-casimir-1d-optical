using QuadGK
using Plots
using LaTeXStrings
using LinearAlgebra

## Parameters
const ν = 1e-2;         # Small number for relative tolerance
const α = 1e-5;         # Small number for absolute tolerance
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
function F_I_Integrand(z, imps)
    nImps = length(imps)
    Λ_ = 1 ./ map(x -> x.M, imps) .- 1 |> Diagonal
    Δ_ = map(x -> x.Δ, imps) |> Diagonal
    Imp_Mat = repeat(imps, 1, nImps)    # Impurity matrix
    ImpT_Mat = permutedims(Imp_Mat)     # Transpose of impurity matrix
    Y = map((x, y) -> propagator(z, abs.(x.pos - y.pos)), Imp_Mat, ImpT_Mat)
    return (
               Matrix{Int}(I, 2 * nImps, 2 * nImps) + [
                   Y * Δ_ z .* Y * Λ_
                   z .* Y * Δ_ (Matrix{Int}(I, nImps, nImps) + z^2 .* Y) * Λ_
               ]
           ) |>
           det |>
           log
end

function F_I(imps)
    f(x) = F_I_Integrand(1im * x, imps)
    res = quadgk(f, 0, Inf, maxevals = NumEvals)[1]
    return (res / (2 * π))
end

function Exact_Free_Energy(N, imps)
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
    return (sum(Ω) / 2)
end
