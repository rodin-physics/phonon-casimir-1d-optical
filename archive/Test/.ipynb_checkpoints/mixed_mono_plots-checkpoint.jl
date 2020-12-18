include("../src/mixed_library.jl")

Ds = 1:50
Ds_Exact = 2 : 3 : 50   # Distances used for the exact diagonalization
N = 1000                # Chain length for the exact diagonalization
# T in units of Ω0
T = 1e-12;
Ds = Ds_Exact
# Impurity masses in units of m (chain masses)
scale = [0.0001, 0.01, 1., 100., 1000.]


γs = [3, 3/4, 1/4, 1/10, 1/1000]
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



function F_I_int(t, imp_set)
    imp_m = Impurity(0, imp_set[1], 0)
    imp_t = Impurity(t, 1., imp_set[2])
    imp_t1 = Impurity(1, 1., imp_set[2])
    # F1 = real(F_I([imp_m, imp_t1]) - F_I([imp_m]) - F_I([imp_t1]))
    F1 = 1.
    return real(F_I([imp_m, imp_t]) - F_I([imp_m]) - F_I([imp_t])) / F1
end

function F_I_int_norm(t, imp_set)
    imp_m = Impurity(0, imp_set[1], 0)
    imp_t = Impurity(t, 1., imp_set[2])
    imp_t1 = Impurity(1, 1., imp_set[2])
    F1 = real(F_I([imp_m, imp_t1]) - F_I([imp_m]) - F_I([imp_t1]))
    return real(F_I([imp_m, imp_t]) - F_I([imp_m]) - F_I([imp_t])) / F1
end

set_1a = [1.2, 1.44]
set_1b = [1.2, 100.]
set_1c = [100., 120.]
set_1d = [100., 10000.]

range_1a = map(t -> F_I_int(t, set_1a), Ds)
range_1b = map(t -> F_I_int(t, set_1b), Ds)
range_1c = map(t -> F_I_int(t, set_1c), Ds)
range_1d = map(t -> F_I_int(t, set_1d), Ds)

log_range_1a = map(t -> log(F_I_int_norm(t, set_1a)), Ds)
log_range_1b = map(t -> log(F_I_int_norm(t, set_1b)), Ds)
log_range_1c = map(t -> log(F_I_int_norm(t, set_1c)), Ds)
log_range_1d = map(t -> log(F_I_int_norm(t, set_1d)), Ds)

pyplot()
plot(
    xaxis = (L"$D$", font(14, "Serif")),
    yaxis = (L"$F_I$", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomright,
    xlims = (0, 50),
    ylims = (-0.03, 0)
    )

plot!(Ds, range_1a, label = string(set_1a))
plot!(Ds, range_1b, label = string(set_1b))
plot!(Ds, range_1c, label = string(set_1c))
plot!(Ds, range_1d, label = string(set_1d))

savefig("mixed_pair_set_1.png")

plot(
    xaxis = (L"$\ln\, D$", font(14, "Serif")),
    yaxis = (L"$\ln \left(F_I/F_I^1\right) $", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :topright,
    ylims = (-10, 0)
    )

plot!(log.(Ds), log_range_1a, label = string(set_1a))
plot!(log.(Ds), log_range_1b, label = string(set_1b))
plot!(log.(Ds), log_range_1c, label = string(set_1c))
plot!(log.(Ds), log_range_1d, label = string(set_1d))

savefig("mixed_pair_set_1_log.png")




set_2a = reverse(set_1a)
set_2b = reverse(set_1b)
set_2c = reverse(set_1c)
set_2d = reverse(set_1d)

range_2a = map(t -> F_I_int(t, set_2a), Ds)
range_2b = map(t -> F_I_int(t, set_2b), Ds)
range_2c = map(t -> F_I_int(t, set_2c), Ds)
range_2d = map(t -> F_I_int(t, set_2d), Ds)

pyplot()
plot(
    xaxis = (L"$D$", font(14, "Serif")),
    yaxis = (L"$F_I$", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomright,
    xlims = (0, 50),
    ylims = (-0.03, 0)
    )

plot!(Ds, range_2a, label = string(set_2a))
plot!(Ds, range_2b, label = string(set_2b))
plot!(Ds, range_2c, label = string(set_2c))
plot!(Ds, range_2d, label = string(set_2d))

savefig("mixed_pair_set_2.png")

log_range_2a = map(t -> log(F_I_int_norm(t, set_2a)), Ds)
log_range_2b = map(t -> log(F_I_int_norm(t, set_2b)), Ds)
log_range_2c = map(t -> log(F_I_int_norm(t, set_2c)), Ds)
log_range_2d = map(t -> log(F_I_int_norm(t, set_2d)), Ds)

plot(
    xaxis = (L"$\ln\, D$", font(14, "Serif")),
    yaxis = (L"$\ln \left(F_I/F_I^1\right) $", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    ylims = (-10, 0)
    )

plot!(log.(Ds), log_range_2a, label = string(set_2a))
plot!(log.(Ds), log_range_2b, label = string(set_2b))
plot!(log.(Ds), log_range_2c, label = string(set_2c))
plot!(log.(Ds), log_range_2d, label = string(set_2d))

savefig("mixed_pair_set_2_log.png")



set_3a = [0.8, 1.2]
set_3b = [0.8, 100.]
set_3c = [0.01, 1.2]
set_3d = [0.01, 100.]

range_3a = map(t -> F_I_int(t, set_3a), Ds)
range_3b = map(t -> F_I_int(t, set_3b), Ds)
range_3c = map(t -> F_I_int(t, set_3c), Ds)
range_3d = map(t -> F_I_int(t, set_3d), Ds)

pyplot()
plot(
    xaxis = (L"$D$", font(14, "Serif")),
    yaxis = (L"$F_I$", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :topright,
    xlims = (0, 50),
    ylims = (0, 0.004)
    )

plot!(Ds, range_3a, label = string(set_3a))
plot!(Ds, range_3b, label = string(set_3b))
plot!(Ds, range_3c, label = string(set_3c))
plot!(Ds, range_3d, label = string(set_3d))

savefig("mixed_pair_set_3.png")

log_range_3a = map(t -> log(F_I_int_norm(t, set_3a)), Ds)
log_range_3b = map(t -> log(F_I_int_norm(t, set_3b)), Ds)
log_range_3c = map(t -> log(F_I_int_norm(t, set_3c)), Ds)
log_range_3d = map(t -> log(F_I_int_norm(t, set_3d)), Ds)

plot(
    xaxis = (L"$\ln\, D$", font(14, "Serif")),
    yaxis = (L"$\ln \left(F_I/F_I^1\right) $", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    ylims = (-10, 0)
    )

plot!(log.(Ds), log_range_3a, label = string(set_3a))
plot!(log.(Ds), log_range_3b, label = string(set_3b))
plot!(log.(Ds), log_range_3c, label = string(set_3c))
plot!(log.(Ds), log_range_3d, label = string(set_3d))

savefig("mixed_pair_set_3_log.png")



set_4a = reverse(set_3a)
set_4b = reverse(set_3b)
set_4c = reverse(set_3c)
set_4d = reverse(set_3d)

range_4a = map(t -> F_I_int(t, set_4a), Ds)
range_4b = map(t -> F_I_int(t, set_4b), Ds)
range_4c = map(t -> F_I_int(t, set_4c), Ds)
range_4d = map(t -> F_I_int(t, set_4d), Ds)

pyplot()
plot(
    xaxis = (L"$D$", font(14, "Serif")),
    yaxis = (L"$F_I$", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomright,
    xlims = (0, 50),
    ylims = (-0.017, 0)
    )

plot!(Ds, range_4a, label = string(set_4a))
plot!(Ds, range_4b, label = string(set_4b))
plot!(Ds, range_4c, label = string(set_4c))
plot!(Ds, range_4d, label = string(set_4d))

savefig("mixed_pair_set_4.png")

log_range_4a = map(t -> log(F_I_int_norm(t, set_4a)), Ds)
log_range_4b = map(t -> log(F_I_int_norm(t, set_4b)), Ds)
log_range_4c = map(t -> log(F_I_int_norm(t, set_4c)), Ds)
log_range_4d = map(t -> log(F_I_int_norm(t, set_4d)), Ds)

plot(
    xaxis = (L"$\ln\, D$", font(14, "Serif")),
    yaxis = (L"$\ln \left(F_I/F_I^1\right) $", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    ylims = (-10, 0)
    )

plot!(log.(Ds), log_range_4a, label = string(set_4a))
plot!(log.(Ds), log_range_4b, label = string(set_4b))
plot!(log.(Ds), log_range_4c, label = string(set_4c))
plot!(log.(Ds), log_range_4d, label = string(set_4d))

savefig("mixed_pair_set_4_log.png")



set_5a = [0.64, 0.8]
set_5b = [0.008, 0.01]
set_5c = [0.01, 0.8]
set_5d = [0.0001, 0.01]

range_5a = map(t -> F_I_int(t, set_5a), Ds)
range_5b = map(t -> F_I_int(t, set_5b), Ds)
range_5c = map(t -> F_I_int(t, set_5c), Ds)
range_5d = map(t -> F_I_int(t, set_5d), Ds)

pyplot()
plot(
    xaxis = (L"$D$", font(14, "Serif")),
    yaxis = (L"$F_I$", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :topright,
    xlims = (0, 50),
    ylims = (0, 0.0025)
    )

plot!(Ds, range_5a, label = string(set_5a))
plot!(Ds, range_5b, label = string(set_5b))
plot!(Ds, range_5c, label = string(set_5c))
plot!(Ds, range_5d, label = string(set_5d))

savefig("mixed_pair_set_5.png")

log_range_5a = map(t -> log(F_I_int_norm(t, set_5a)), Ds)
log_range_5b = map(t -> log(F_I_int_norm(t, set_5b)), Ds)
log_range_5c = map(t -> log(F_I_int_norm(t, set_5c)), Ds)
log_range_5d = map(t -> log(F_I_int_norm(t, set_5d)), Ds)

plot(
    xaxis = (L"$\ln\, D$", font(14, "Serif")),
    yaxis = (L"$\ln \left(F_I/F_I^1\right) $", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    ylims = (-10, 0)
    )

plot!(log.(Ds), log_range_5a, label = string(set_5a))
plot!(log.(Ds), log_range_5b, label = string(set_5b))
plot!(log.(Ds), log_range_5c, label = string(set_5c))
plot!(log.(Ds), log_range_5d, label = string(set_5d))

savefig("mixed_pair_set_5_log.png")



set_6a = [0.8, 1.2]
set_6b = [0.8, 100.]
set_6c = [0.01, 1.2]
set_6d = [0.01, 100.]

range_6a = map(t -> F_I_int(t, set_6a), Ds)
range_6b = map(t -> F_I_int(t, set_6b), Ds)
range_6c = map(t -> F_I_int(t, set_6c), Ds)
range_6d = map(t -> F_I_int(t, set_6d), Ds)

pyplot()
plot(
    xaxis = (L"$D$", font(14, "Serif")),
    yaxis = (L"$F_I$", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :topright,
    xlims = (0, 50),
    ylims = (0, 0.003)
    )

plot!(Ds, range_6a, label = string(set_6a))
plot!(Ds, range_6b, label = string(set_6b))
plot!(Ds, range_6c, label = string(set_6c))
plot!(Ds, range_6d, label = string(set_6d))

savefig("mixed_pair_set_6.png")

log_range_6a = map(t -> log(F_I_int_norm(t, set_6a)), Ds)
log_range_6b = map(t -> log(F_I_int_norm(t, set_6b)), Ds)
log_range_6c = map(t -> log(F_I_int_norm(t, set_6c)), Ds)
log_range_6d = map(t -> log(F_I_int_norm(t, set_6d)), Ds)

plot(
    xaxis = (L"$\ln\, D$", font(14, "Serif")),
    yaxis = (L"$\ln \left(F_I/F_I^1\right) $", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    ylims = (-10, 0)
    )

plot!(log.(Ds), log_range_6a, label = string(set_6a))
plot!(log.(Ds), log_range_6b, label = string(set_6b))
plot!(log.(Ds), log_range_6c, label = string(set_6c))
plot!(log.(Ds), log_range_6d, label = string(set_6d))

savefig("mixed_pair_set_6_log.png")






Ω_0 = 1.
γ = 1 / 1.2
λ = 1 / 1.44
t = 16

# seems like impurity mass should be 1
test_imps = [Impurity(0, 1 / γ, η), Impurity(t, η, λ)]
imp_m = Impurity(0, 1 / γ, 0)
imp_t(t) = Impurity(t, 1, λ)

F_I_int(t, λ, 1 / γ)

F_I([imp_m, imp_t(t)]) - F_I([imp_m, Impurity(t, 1 + η, λ)])
F_I([imp_m, imp_t(t)]) - F_I([imp_m]) - F_I([imp_t(t)])
F_I([imp_t(t)]) - F_I([imp_t(t), Impurity(0, 1, 0)])


F_I([imp_m, imp_t])
F_0 = F_I([Impurity(0, η, η)])

range_of_F = map(t -> real(F_I([imp_m, imp_t(t)]) - F_I([imp_t(t)]) - F_I([imp_m])), Ds)
plot(Ds, range_of_F)

F_0 = Exact_Free_Energy(2000, [])
@time Exact_Free_Energy(2000, [imp_m, imp_t(t)]) + F_0 - Exact_Free_Energy(2000, [imp_t(t)]) - Exact_Free_Energy(2000, [imp_m])

F_0_1000 = Exact_Free_Energy(1000, [])
@time Exact_Free_Energy(1000, [imp_m, imp_t(t)]) + F_0_1000 - Exact_Free_Energy(1000, [imp_t(t)]) - Exact_Free_Energy(1000, [imp_m])

range_of_F = map(t -> real(F_I([Impurity(0, 1 / γ, η), Impurity(t, η, λ)]) - F_I([Impurity(0, 1 / γ, η)]) - F_I([Impurity(t, η, λ)])), Ds)
plot(Ds, range_of_F)

Exact_Free_Energy(2000, [Impurity(0, 1 / γ, η), Impurity(t, η, λ)]) - Exact_Free_Energy(2000, [Impurity(0, 1 / γ, η)]) - Exact_Free_Energy(2000, [Impurity(0, 1 / γ, η)])
@time range_of_F_exact = map(t -> Exact_Free_Energy(1000, [imp_m, imp_t(t)]) + F_0_1000 - Exact_Free_Energy(1000, [imp_t(t)]) - Exact_Free_Energy(1000, [imp_m]), Ds_Exact)
plot!(Ds_Exact, range_of_F_exact)

F_I([imp_m, imp_t])


Exact_Free_Energy(2000, [Impurity(0, 1 / γ, 0), Impurity(t, η, λ)])  - Exact_Free_Energy(2000, [Impurity(0, 1 / γ, 0), Impurity(t, η * 0.1, λ)])

N = 1000

E_0 = Exact_Free_Energy(N, [])
E_I = Exact_Free_Energy(N, [imp_m, imp_t]) - E_0
E_m = Exact_Free_Energy(N, [imp_m]) - E_0
E_t = Exact_Free_Energy(N, [imp_t]) - E_0
E_int = E_I - E_m - E_t
x = Exact_Free_Energy(2000, [imp_m, imp_t]) - Exact_Free_Energy(2000, [imp_m]) - Exact_Free_Energy(2000, [imp_t])
x




m = Ms[1]
Ω_0 = 1.
γ = γs[1]
λ = λs[1]
t = 16

function F_I_int2(d, imp_set)
    imp_t_L = Impurity(0, 1., imp_set[1])
    imp_t_R = Impurity(30, 1., imp_set[2])
    imp_m = Impurity(d, imp_set[3], 0)
    return real(F_I([imp_t_L, imp_t_R, imp_m]) - F_I([imp_m]) - F_I([imp_t_L]) -F_I([imp_t_R]))
end

function F_I_int2_norm(d, imp_set)
    imp_t_L = Impurity(0, 1., imp_set[1])
    imp_t_R = Impurity(30, 1., imp_set[2])
    imp_m = Impurity(d, imp_set[3], 0)
    imp_m0 = Impurity(15, imp_set[3], 0)
    Fmid = real(F_I([imp_t_L, imp_t_R, imp_m0]) - F_I([imp_m0]) - F_I([imp_t_L]) -F_I([imp_t_R]))
    Freal = real(F_I([imp_t_L, imp_t_R, imp_m]) - F_I([imp_m]) - F_I([imp_t_L]) -F_I([imp_t_R]))
    return  abs((Fmid - Freal) / Fmid)
end

Ds_2 = 1:29
half_Ds_2 = 1:14

set2_1a = [15., 15., 3.]
set2_1b = [15., 15., 0.8]

range2_1a = map(d -> F_I_int2(d, set2_1a), Ds_2)
range2_1b = map(d -> F_I_int2(d, set2_1b), Ds_2)
plot(Ds_2, range2_1a, label = string(set2_1a))
plot!(Ds_2, range2_1b, label = string(set2_1b))

savefig("mixed_triple_set_1.png")

log_range2_1a = map(d -> log(F_I_int2_norm(d, set2_1a)), Ds_2)
log_range2_1b = map(d -> log(F_I_int2_norm(d, set2_1b)), Ds_2)
plot(log.(half_Ds_2), log_range2_1a[16:29], label = string(set2_1a, "_R"))
plot!(log.(half_Ds_2), reverse(log_range2_1a[1:14]), label = string(set2_1a, "_L"))
plot!(log.(half_Ds_2), log_range2_1b[16:29], label = string(set2_1b, "_R"))
plot!(log.(half_Ds_2), reverse(log_range2_1b[1:14]), label = string(set2_1b, "_L"))

savefig("mixed_triple_set_1_log.png")

set2_2a = [15., 5., 3.]
set2_2b = [15., 5., 0.8]

range2_2a = map(d -> F_I_int2(d, set2_2a), Ds_2)
range2_2b = map(d -> F_I_int2(d, set2_2b), Ds_2)
plot(Ds_2, range2_2a, label = string(set2_2a))
plot!(Ds_2, range2_2b, label = string(set2_2b))

savefig("mixed_triple_set_2.png")

log_range2_2a = map(d -> log(F_I_int2_norm(d, set2_2a)), Ds_2)
log_range2_2b = map(d -> log(F_I_int2_norm(d, set2_2b)), Ds_2)
plot(log.(half_Ds_2), log_range2_2a[16:29], label = string(set2_2a, "_R"))
plot!(log.(half_Ds_2), reverse(log_range2_2a[1:14]), label = string(set2_2a, "_L"))
plot!(log.(half_Ds_2), log_range2_2b[16:29], label = string(set2_2b, "_R"))
plot!(log.(half_Ds_2), reverse(log_range2_2b[1:14]), label = string(set2_2b, "_L"))

savefig("mixed_triple_set_2_log.png")

set2_3a = [100000., 5., 3.]
set2_3b = [100000., 5., 0.8]

range2_3a = map(d -> F_I_int2(d, set2_3a), Ds_2)
range2_3b = map(d -> F_I_int2(d, set2_3b), Ds_2)
plot(Ds_2, range2_3a, label = string(set2_3a))
plot!(Ds_2, range2_3b, label = string(set2_3b))

savefig("mixed_triple_set_3.png")

log_range2_3a = map(d -> log(F_I_int2_norm(d, set2_3a)), Ds_2)
log_range2_3b = map(d -> log(F_I_int2_norm(d, set2_3b)), Ds_2)
plot(log.(half_Ds_2), log_range2_3a[16:29], label = string(set2_3a, "_R"))
plot!(log.(half_Ds_2), reverse(log_range2_3a[1:14]), label = string(set2_3a, "_L"))
plot!(log.(half_Ds_2), log_range2_3b[16:29], label = string(set2_3b, "_R"))
plot!(log.(half_Ds_2), reverse(log_range2_3b[1:14]), label = string(set2_3b, "_L"))

savefig("mixed_triple_set_3_log.png")


# set2_4a = [15., 1.5, 3.]
# set2_4b = [15., 0.5, 0.8]
set2_4c = [15., 0.0001, 3.]
set2_4d = [15., 0.0001, 0.8]

# range2_4a = map(d -> F_I_int2(d, set2_4a), Ds_2)
# range2_4b = map(d -> F_I_int2(d, set2_4b), Ds_2)
range2_4c = map(d -> F_I_int2(d, set2_4c), Ds_2)
range2_4d = map(d -> F_I_int2(d, set2_4d), Ds_2)
# plot(Ds_2, range2_4a, label = string(set2_4a))
# plot!(Ds_2, range2_4b, label = string(set2_4b))
plot(Ds_2, range2_4c, label = string(set2_4c))
plot!(Ds_2, range2_4d, label = string(set2_4d))

savefig("mixed_triple_set_4.png")

log_range2_4c = map(d -> log(F_I_int2_norm(d, set2_4c)), Ds_2)
log_range2_4d = map(d -> log(F_I_int2_norm(d, set2_4d)), Ds_2)
plot(log.(half_Ds_2), log_range2_4c[16:29], label = string(set2_4c, "_R"))
plot!(log.(half_Ds_2), reverse(log_range2_4c[1:14]), label = string(set2_4c, "_L"))
plot!(log.(half_Ds_2), log_range2_4d[16:29], label = string(set2_4d, "_R"))
plot!(log.(half_Ds_2), reverse(log_range2_4d[1:14]), label = string(set2_4d, "_L"))

savefig("mixed_triple_set_4_log.png")

λ
