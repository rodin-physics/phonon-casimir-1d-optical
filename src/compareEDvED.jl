include("monoatomic_library.jl")
include("exactDiag1D.jl")

# Comparing for the 0 temperature case

Ds = 2:50
Ds_Exact = 2 : 3 : 50   # Distances used for the exact diagonalization
N = 1000                # Chain length for the exact diagonalization
# T in units of Ω0
T = 1e-12;

# Impurity masses in units of m (chain masses)
Ms = [1/3, 4/3, 4, 10, 1000]

# Bounding power laws
r_1 = 1 ./ Ds
r_3 = 1 ./ (Ds.^3)

# special case definition so the two kinds of impurities don't clash
struct diImpurity
    pos::Int    # Position of the impurity unit cell
    n::Int      # Index of the atom in the unit cell
    λ::Float64  # Impurity mass in units of μ
end

## Plotting
pyplot()
plot(
    xaxis = (L"$\ln\, D$", font(14, "Serif")),
    yaxis = (L"$\ln \left(F_I/F_I^1\right) $", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    yticks = -10 : 5 : 0
    )

## Baseline
plot!(
    log.(Ds), log.(r_1),
    linewidth = 2,
    color = RGB(50/255,50/255,50/255),
    line = :dash,
    lab = ""
    )
plot!(
    log.(Ds), log.(r_3),
    linewidth = 2,
    color = RGB(50/255,50/255,50/255),
    line = :dash,
    lab = ""
    )

## Reference plot using monoatomic_library
for ii = 1 : (length(Ms) - 0)
    println(ii)
    M = Ms[ii];
    # Energy for maximally separated impurities in a finite-length chain
    E_halfway = Exact_Free_Energy(N, M, M, floor(Int, N / 2), T)
    # F_I for adjacent impurities
    E0 = Exact_Free_Energy(N, M, M, 1, T) - E_halfway
    # F_I divided by F_I at D = 1
    r =  map(x -> Exact_Free_Energy(N, M, M, x, T) - E_halfway, Ds_Exact) ./ E0;
    Plots.display(plot!(log.(Ds_Exact), real(log.(complex(r))),
        color = colors[ii],
        lab = "",
        markershape = :circle
        ))
end

function two_diImps(D, impM)
    Imps = [diImpurity(1, 1, impM), diImpurity(1 + floor(Int, D/2), 1 + Int(D % 2), impM)]
    return Imps
end

for ii = 1: (length(Ms) - 0)
    diM = [1., 1.];
    M = Ms[ii];
    F_halfway = exact_F(diM, floor(Int, N / 2), two_diImps(floor(Int, N/2), M), T)
    F0 = exact_F(diM, floor(Int, N / 2), two_diImps(1, M), T) - F_halfway
    r =  map(x -> exact_F(diM, floor(Int, N / 2), two_diImps(x, M), T) - F_halfway, Ds_Exact) ./ F0;

    Plots.display(plot!(log.(Ds_Exact), real(log.(complex(r))),
        color = colors[ii],
        lab = "",
        markershape = :cross
        ))

end

## Test plot using exactDiag1D
