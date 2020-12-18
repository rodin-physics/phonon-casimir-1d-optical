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

function dprop_finitediff(z, D, dz)
    println(propagator(z + dz, D))
    println(propagator(z, D))
    return (propagator(z + dz, D) - propagator(z, D)) / dz
end

function dprop_analytic(z, D)
    # x1 = 1 - 2 * z^2 - 2 * sqrt(z ^ 4 - z ^ 2)
    # x2 = (4 * z^3 - 2 * z) / sqrt(z ^ 4 - z ^ 2)
    if D > 0
        x1 = 1 - 2 * z^2 - 2 * sqrt( - z^2) * sqrt(1 - z^2)
        x2 = sqrt( - z^2) * ( 1 - 2 * z^2) * sqrt(1 - z^2) + 2 * D * z^2 * (-1 + z^2)
        x3 = z^3 * (-1 + z^2)^2
        upper = x1 ^ D * x2
        lower = x3
        # upper = D * x1 ^ (D - 1) * (- 4 * z - x2) - (1 / 2) * x1 ^ D * x2
        # lower = z^4 - z^2
    elseif D == 0
        upper = z * (1 - 2 * z ^ 2)
        lower = (-z^2) ^ (3 / 2) * (1 - z^2 ) ^ (3/2)
    end
    return upper / lower
end

# checking differentiating Y at D = 0 ... works pretty well!
D = 0
z = 0 - 4im
z = 4 + 0im
z = 2 + 2im
dz = 0 + 1e-10im
dprop_analytic(z, D)
dprop_finitediff(z, D, dz)
dz = 1e-10 + 0im
dprop_analytic(z, D)
dprop_finitediff(z, D, dz)
dz = 1e-10 + 1e-10im
dprop_analytic(z, D)
dprop_finitediff(z, D, dz)

# the mathematica formula works well for finite D as well
D = 17
z = 0 - 4im
z = 4 + 0im
z = 2 + 2im
dz = 0 + 1e-10im
dprop_analytic(z, D)
dprop_finitediff(z, D, dz)
dz = 1e-10 + 0im
dprop_analytic(z, D)
dprop_finitediff(z, D, dz)
dz = 1e-10 + 1e-10im
dprop_analytic(z, D)
dprop_finitediff(z, D, dz)

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


function dAz_finitediff(z, imps, dz)
    return (Az(z + dz, imps) - Az(z, imps)) ./ dz
end

function dAz_analytic(z, imps)
    nImps = length(imps)
    Λ_ = 1 ./ map(x -> x.M, imps) .- 1 |> Diagonal
    Δ_ = map(x -> x.Δ, imps) |> Diagonal
    Imp_Mat = repeat(imps, 1, nImps)    # Impurity matrix
    ImpT_Mat = permutedims(Imp_Mat)     # Transpose of impurity matrix
    Y = map((x, y) -> propagator(z, abs.(x.pos - y.pos)), Imp_Mat, ImpT_Mat)
    dY = map((x, y) -> dprop_analytic(z, abs.(x.pos - y.pos)), Imp_Mat, ImpT_Mat)
    return [
        dY * Δ_ (Y + z .* dY) * Λ_
        (Y + z .* dY) * Δ_ (2 * z .* Y + z^2 .* dY) * Λ_
    ]
end

imp_m = Impurity(0, 3., 0)
imp_t = Impurity(17, 1., 1/3)
testimps = [imp_m, imp_t]
Az(0.5 + 1e-7im, testimps)

# dZ seems to work pretty well too!
z = 0 - 4im
z = 4 + 0im
z = 2 + 2im
dz = 0 + 1e-10im
dAz_analytic(z, testimps)[1, 1]
dAz_finitediff(z, testimps, dz)[1, 1]
dAz_analytic(z, testimps)[1, 2]
dAz_finitediff(z, testimps, dz)[1, 2]
dAz_analytic(z, testimps)[2, 1]
dAz_finitediff(z, testimps, dz)[2, 1]
dAz_analytic(z, testimps)[2, 2]
dAz_finitediff(z, testimps, dz)[2, 2]
dz = 1e-10 + 0im
dAz_analytic(z, testimps)[1, 1]
dAz_finitediff(z, testimps, dz)[1, 1]
dAz_analytic(z, testimps)[1, 2]
dAz_finitediff(z, testimps, dz)[1, 2]
dAz_analytic(z, testimps)[2, 1]
dAz_finitediff(z, testimps, dz)[2, 1]
dAz_analytic(z, testimps)[2, 2]
dAz_finitediff(z, testimps, dz)[2, 2]
dz = 1e-10 + 1e-10im
dAz_analytic(z, testimps)[1, 1]
dAz_finitediff(z, testimps, dz)[1, 1]
dAz_analytic(z, testimps)[1, 2]
dAz_finitediff(z, testimps, dz)[1, 2]
dAz_analytic(z, testimps)[2, 1]
dAz_finitediff(z, testimps, dz)[2, 1]
dAz_analytic(z, testimps)[2, 2]
dAz_finitediff(z, testimps, dz)[2, 2]

function F_I_Integrand(z, imps, T)
    nImps = length(imps)
    Az_mat = Az(z, imps)
    dAz_mat = dAz_analytic(z, imps)
    term1 = log(exp(- z / T))
    term2 = tr(inv(Matrix{Int}(I, 2 * nImps, 2 * nImps) + Az_mat) * dAz_mat)
    return term1 * term2
end

function naive_Summand(z, imps, T)
    nImps = length(imps)
    Az_mat = Az(z, imps)
    return (T / 2) * log(det(1 + Az_mat))
end

F_I_Integrand(0.5, testimps, T)

T = 5. + 0im
F_I_Integrand(z, testimps, T)
coord_grid = hcat(-2:2, -2:2, -2:2, -2:2, -2:2)
test_vals = map((x, y) -> abs.(F_I_Integrand(x + y * 1im, testimps, T)), transpose(coord_grid), coord_grid)
gr()
print(typeof(test_vals))
print(test_vals)
print(size(test_vals))
heatmap(1:size(test_vals,1),
    1:size(test_vals,2), test_vals,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="x values", ylabel="y values",
    title="My title")

heatmap(-5:0.5:5,
    -5:0.5:5, test_vals,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="x values", ylabel="y values",
    title="My title")

heatmap(-5:5:10,
    -5:5:10, test_vals,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="x values", ylabel="y values",
    title="My title")


for ii = -5:0.5:5
    print(ii)
end

using Plots
gr()
data = rand(21,100)
heatmap(1:size(data,1),
    1:size(data,2), data,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="x values", ylabel="y values",
    title="My title")


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
