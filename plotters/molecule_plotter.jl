include("../calculations/molecule_calc_settings.jl")
pgfplotsx()
## Zero-T plots

# F_I as a function of M
exact_data = readdlm("data/molecule/F_I_exact_T0_M.dat")
QFT_data = readdlm("data/molecule/F_I_QFT_T0_M.dat")

x_min = M_min
x_max = M_max
y_min = -0.5
y_max = 0.025

plot(
    xaxis = (L"M / m", font(14, "Serif")),
    yaxis = (L"F_I/\Omega", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    ylims = (y_min, y_max),
    xlims = (x_min, x_max),
    legend = :bottomright,
)

for ii = 1:length(Ks)
    # Plot the ExD data
    plot!(
        Ms_exact,
        exact_data[:, ii],
        linewidth = 2,
        color = colors[ii],
        label = false,
    )
    # Plot the QFT data
    scatter!(
        Ms_QFT,
        QFT_data[:, ii],
        label = latexstring("\$K/k = $(Ks[ii])\$"),
        color = colors[ii],
    )
    # Plot the asymptote
    Ω0_a = √(Ks[ii])
    Ω0_b = √(2 + Ks[ii])
    Ω1 = √(1 + Ks[ii])
    F0 = (Ω0_a + Ω0_b) / 2
    F1 = Ω1 / 2
    F_I_asymptote = F0 - 2 * F1

    plot!(
        [M_min, M_max],
        F_I_asymptote .* [1, 1],
        linestyle = :dot,
        linewidth = 2,
        lab = false,
        color = colors[ii],
    )
end
# K = 0 asymptote
plot!(
    [M_min, M_max],
    (√(2) - 2) / 2 .* [1, 1],
    linestyle = :dash,
    linewidth = 2,
    lab = false,
    color = RGB(50 / 255, 50 / 255, 50 / 255),
)
annotate!(
    x_min + (0.1) * (x_max - x_min),
    y_min + (0.075) * (y_max - y_min),
    text("(a)", 18),
)

savefig("plotters/plots/molecule/Dimer_T0_M.pdf")

# F_I as a function of ln(M)
exact_data = readdlm("data/molecule/F_I_exact_T0_log_M.dat")
QFT_data = readdlm("data/molecule/F_I_QFT_T0_log_M.dat")

x_min = log(M_min)
x_max = log(M_max)
y_min = -0.5
y_max = 0.025

plot(
    xaxis = (L"\ln(M / m)", font(14, "Serif")),
    yaxis = (L"F_I/\Omega", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    ylims = (y_min, y_max),
    xlims = (x_min, x_max),
    legend = :bottomright,
)

for ii = 1:length(Ks)
    # Plot the ExD data
    plot!(
        log.(log_Ms_exact),
        exact_data[:, ii],
        linewidth = 2,
        color = colors[ii],
        label = false,
    )
    # Plot the QFT data
    scatter!(
        log.(log_Ms_QFT),
        QFT_data[:, ii],
        label = latexstring("\$K/k = $(Ks[ii])\$"),
        color = colors[ii],
    )
    # Plot the asymptote
    Ω0_a = √(Ks[ii])
    Ω0_b = √(2 + Ks[ii])
    Ω1 = √(1 + Ks[ii])
    F0 = (Ω0_a + Ω0_b) / 2
    F1 = Ω1 / 2
    F_I_asymptote = F0 - 2 * F1

    plot!(
        log.([M_min, M_max]),
        F_I_asymptote .* [1, 1],
        linestyle = :dot,
        linewidth = 2,
        lab = false,
        color = colors[ii],
    )
end
# K = 0 asymptote
plot!(
    log.([M_min, M_max]),
    (√(2) - 2) / 2 .* [1, 1],
    linestyle = :dash,
    linewidth = 2,
    lab = false,
    color = RGB(50 / 255, 50 / 255, 50 / 255),
)
annotate!(
    x_min + (0.1) * (x_max - x_min),
    y_min + (0.075) * (y_max - y_min),
    text("(b)", 18),
)

savefig("plotters/plots/molecule/Dimer_T0_Log_M.pdf")


# F_I as a function of δ
exact_data = readdlm("data/molecule/F_I_exact_T0_delta.dat")
QFT_data = readdlm("data/molecule/F_I_QFT_T0_delta.dat")

x_min = δ_min
x_max = δ_max
y_min = -0.5
y_max = 0.025

plot(
    xaxis = (L"\delta / k", font(14, "Serif")),
    yaxis = (L"F_I/\Omega", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    ylims = (y_min, y_max),
    xlims = (x_min, x_max),
    legend = :bottomright,
)

for ii = 1:length(Ks)
    # Plot the ExD data
    plot!(
        δs_exact,
        exact_data[:, ii],
        linewidth = 2,
        color = colors[ii],
        label = false,
    )
    # Plot the QFT data
    scatter!(
        δs_QFT,
        QFT_data[:, ii],
        label = latexstring("\$K/k = $(Ks[ii])\$"),
        color = colors[ii],
    )
    # Plot the asymptote
    Ω0_a = √(Ks[ii])
    Ω0_b = √(2 + Ks[ii])
    Ω1 = √(1 + Ks[ii])
    F0 = (Ω0_a + Ω0_b) / 2
    F1 = Ω1 / 2
    F_I_asymptote = F0 - 2 * F1

    plot!(
        [δ_min, δ_max],
        F_I_asymptote .* [1, 1],
        linestyle = :dot,
        linewidth = 2,
        lab = false,
        color = colors[ii],
    )
end
# K = 0 asymptote
plot!(
    [δ_min, δ_max],
    (√(2) - 2) / 2 .* [1, 1],
    linestyle = :dash,
    linewidth = 2,
    lab = false,
    color = RGB(50 / 255, 50 / 255, 50 / 255),
)
annotate!(
    x_min + (0.1) * (x_max - x_min),
    y_min + (0.075) * (y_max - y_min),
    text("(c)", 18),
)

savefig("plotters/plots/molecule/Dimer_T0_Delta.pdf")

## Finite-T plots

# F_I fixed M
exact_data = readdlm("data/molecule/F_I_exact_T_single_M.dat")
QFT_data = readdlm("data/molecule/F_I_QFT_T_single_M.dat")

x_min = T_min - 0.05
x_max = T_max + 0.05
y_min = -0.125
y_max = 0.005

plot(
    xaxis = (L"T / \Omega", font(14, "Serif")),
    yaxis = (L"F_I/\Omega", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    ylims = (y_min, y_max),
    xlims = (x_min, x_max),
    legend = :bottomright,
    yticks = [-0.1, 0],

)

for ii = 1:length(Ks)
    # Plot the ExD data
    plot!(
        Ts_exact,
        exact_data[:, ii],
        linewidth = 2,
        color = colors[ii],
        label = false,
    )
    # Plot the QFT data
    scatter!(
        Ts_QFT,
        QFT_data[:, ii],
        label = latexstring("\$K/k = $(Ks[ii])\$"),
        color = colors[ii],
    )
end

annotate!(
    x_min + (0.1) * (x_max - x_min),
    y_min + (0.075) * (y_max - y_min),
    text("(a)", 18),
)

savefig("plotters/plots/molecule/Dimer_T_single_M.pdf")

# F_I fixed δ
exact_data = readdlm("data/molecule/F_I_exact_T_single_delta.dat")
QFT_data = readdlm("data/molecule/F_I_QFT_T_single_delta.dat")

x_min = T_min - 0.05
x_max = T_max + 0.05
y_min = -2.5
y_max = 0.005

plot(
    xaxis = (L"T / \Omega", font(14, "Serif")),
    yaxis = (L"F_I/\Omega", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    ylims = (y_min, y_max),
    xlims = (x_min, x_max),
    legend = :bottomleft,
    # yticks = [-0.1, 0],

)

for ii = 1:length(Ks)
    # Plot the ExD data
    plot!(
        Ts_exact,
        exact_data[:, ii],
        linewidth = 2,
        color = colors[ii],
        label = false,
    )
    # Plot the QFT data
    scatter!(
        Ts_QFT,
        QFT_data[:, ii],
        label = latexstring("\$K/k = $(Ks[ii])\$"),
        color = colors[ii],
    )
end

annotate!(
    x_min + (0.95) * (x_max - x_min),
    y_min + (0.075) * (y_max - y_min),
    text("(b)", 18),
)

savefig("plotters/plots/molecule/Dimer_T_single_Delta.pdf")

# F_I fixed δ
exact_data = readdlm("data/molecule/E_I_exact_T_single_delta.dat")
QFT_data = readdlm("data/molecule/E_I_QFT_T_single_delta.dat")

x_min = T_min - 0.05
x_max = T_max + 0.05
y_min = -.125
y_max = 0.005

plot(
    xaxis = (L"T / \Omega", font(14, "Serif")),
    yaxis = (L"(F_I+TS_I)/\Omega", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    ylims = (y_min, y_max),
    xlims = (x_min, x_max),
    legend = :bottomright,
    yticks = [-0.1, 0],

)

for ii = 1:length(Ks)
    # Plot the ExD data
    plot!(
        Ts_exact,
        exact_data[:, ii],
        linewidth = 2,
        color = colors[ii],
        label = false,
    )
    # Plot the QFT data
    scatter!(
        Ts_QFT,
        QFT_data[:, ii],
        label = latexstring("\$K/k = $(Ks[ii])\$"),
        color = colors[ii],
    )
end

annotate!(
    x_min + (0.1) * (x_max - x_min),
    y_min + (0.075) * (y_max - y_min),
    text("(c)", 18),
)

savefig("plotters/plots/molecule/Dimer_T_single_Delta_EI.pdf")
