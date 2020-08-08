include("../src/mixed_library.jl")

Ds = 2:50
Ds_Exact = 2 : 3 : 50   # Distances used for the exact diagonalization
N = 1000                # Chain length for the exact diagonalization
# T in units of Ω0
T = 1e-12;
Ds = Ds_Exact
# Impurity masses in units of m (chain masses)
γs = [1/3, 4/3, 4, 10, 1000]
λs = [1/3, 4/3, 4, 10, 1000]
# Chain masses
Ms = [1/3, 1, 3, 12, 50]

# Bounding power laws
r_1 = 1 ./ Ds
r_3 = 1 ./ (Ds.^3)

pyplot()
plot(
    xaxis = (L"$\ln\, D$", font(14, "Serif")),
    yaxis = (L"$\ln \left(F_I/F_I^1\right) $", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    yticks = -10 : 5 : 0,
    ylims = (-10, 5)
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


m = Ms[1]
Ω_0 = 1.
γ = γs[1]
λ = λs[1]
t = 16

test_imps = [Impurity(0, 1 / γ, η), Impurity(t, η, λ)]
imp_m = Impurity(0, 1 / γ, η)
imp_t = Impurity(t, η, λ)
F_I([imp_m, imp_t]) - F_I([imp_m]) - F_I([imp_t]) + F_0

F_I([imp_m, imp_t])
F_0 = F_I([Impurity(0, η, η)])

range_of_F = map(t -> real(F_I([Impurity(0, 1 / γ, η), Impurity(t, η, λ)]) - F_I([Impurity(0, 1 / γ, η)]) - F_I([Impurity(t, η, λ)])), Ds)
plot(Ds, range_of_F)

Exact_Free_Energy(2000, [Impurity(0, 1 / γ, η), Impurity(t, η, λ)]) - Exact_Free_Energy(2000, [Impurity(0, 1 / γ, η)]) - Exact_Free_Energy(2000, [Impurity(0, 1 / γ, η)])
range_of_F_exact = map(t -> Exact_Free_Energy(1000, [Impurity(0, 1 / γ, η), Impurity(t, η, λ)]) - Exact_Free_Energy(1000, [Impurity(0, 1 / γ, η)]) - Exact_Free_Energy(1000, [Impurity(0, η, )]), Ds_Exact)

F_I([imp_m, imp_t])

N = 1000

E_0 = Exact_Free_Energy(N, [])
E_I = Exact_Free_Energy(N, [imp_m, imp_t]) - E_0
E_m = Exact_Free_Energy(N, [imp_m]) - E_0
E_t = Exact_Free_Energy(N, [imp_t]) - E_0
E_int = E_I - E_m - E_t
x = Exact_Free_Energy(2000, [imp_m, imp_t]) - Exact_Free_Energy(2000, [imp_m]) - Exact_Free_Energy(2000, [imp_t])
x
