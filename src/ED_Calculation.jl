include("general.jl")

# Comparing for the 0 temperature case
Ds = 1:1:50         # Distances used for the analytic calculation
Ds_Exact = 2:3:50   # Distances used for the exact diagonalization
N = 1000          # Chain length for the exact diagonalization
# T in units of Ω0
T = 1e-12;

# Impurity masses in units of m (chain masses)
Ms = [1 / 3, 4 / 3, 4, 10, 1000]
m = 1.0

# Bounding power laws
r_1 = 1 ./ Ds
r_3 = 1 ./ (Ds .^ 3)

## Plotting
pyplot()
plot(
    xaxis = (L"$\ln\, D$", font(14, "Serif")),
    yaxis = (L"$\ln \left(F_I/F_I^1\right) $", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    yticks = -10:5:0,
)

## Baseline
plot!(
    log.(Ds),
    log.(r_1),
    linewidth = 2,
    color = RGB(50 / 255, 50 / 255, 50 / 255),
    line = :dash,
    lab = "",
)

plot!(
    log.(Ds),
    log.(r_3),
    linewidth = 2,
    color = RGB(50 / 255, 50 / 255, 50 / 255),
    line = :dash,
    lab = "",
)

## Use the new exact diagonalization for a monoatomic chain
for ii = 1:length(Ms)
    println(ii)
    M = Ms[ii]

    # Energy for maximally separated impurities in a finite-length chain
    E_halfway = exact_F(
        [m],
        N,
        [Impurity(1, 1, M), Impurity(floor(Int, N / 2), 1, M)],
        T,
    )
    # F_I for adjacent impurities
    E0 = exact_F([m], N, [Impurity(1, 1, M), Impurity(2, 1, M)], T) - E_halfway
    # F_I divided by F_I at D = 1
    r =
        map(
            x ->
                exact_F([m], N, [Impurity(1, 1, M), Impurity(1 + x, 1, M)], T) -
                E_halfway,
            Ds_Exact,
        ) ./ E0
    plot!(
        log.(Ds_Exact),
        log.(r),
        color = colors[ii],
        lab = "",
        markershape = :circle,
    )
end

## Repeat for a diatomic chain with equal masses
for ii = 1:length(Ms)
    println(ii)
    M = Ms[ii]
    # Energy for maximally separated impurities in a finite-length chain
    E_halfway = exact_F(
        [m, m],
        floor(Int, N / 2),
        [Impurity(1, 1, M), Impurity(floor(Int, N / 4), 1, M)],
        T,
    )
    # F_I for adjacent impurities
    E0 =
        exact_F(
            [m, m],
            floor(Int, N / 2),
            [Impurity(1, 1, M), Impurity(1, 2, M)],
            T,
        ) - E_halfway
    # F_I divided by F_I at D = 1
    r =
        map(
            x ->
                exact_F(
                    [m, m],
                    floor(Int, N / 2),
                    [
                        Impurity(1, 1, M),
                        Impurity(1 + floor(Int, x / 2), 1 + isodd(x), M),
                    ],
                    T,
                ) - E_halfway,
            Ds_Exact,
        ) ./ E0
    plot!(
        log.(Ds_Exact),
        log.(r),
        color = colors[ii],
        lab = "",
        markershape = :diamond,
    )
end

savefig("Test.pdf")
