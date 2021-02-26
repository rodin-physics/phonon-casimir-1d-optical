include("../calculations/monoatomic_calc_settings.jl")
pgfplotsx()

## Mass-mass interaction plots
x_min = -0.1;
x_max = log(maximum(ds)) + 0.1;
y_min = -12.5
y_max = 0.25

# T = 0
ExD = readdlm("data/monoatomic_chain/ExD_M-M_T0.dat")
QFT = readdlm("data/monoatomic_chain/QFT_M-M_T0.dat")

plot(
    xaxis = (L"\ln d", font(14, "Serif")),
    yaxis = (L"\ln \left(E_I/E_I^1\right)", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    yticks = -10:5:0,
    xlims = (x_min, x_max),
    ylims = (y_min, y_max),
)
for ii = 1:length(Defect_Ms)
    plot!(
        log.(ds),
        log.(ExD[ii, :] ./ ExD[ii, 1]),
        color = colors[ii],
        label = false,
        linewidth = 2,
    )
    scatter!(
        log.(ds),
        log.(QFT[ii, :] ./ QFT[ii, 1]),
        color = colors[ii],
        label = latexstring("\$M = $(Defect_Ms[ii])\$"),
    )
end

plot!(
    log.(ds),
    log.(1 ./ ds),
    linewidth = 2,
    color = RGB(50 / 255, 50 / 255, 50 / 255),
    line = :dash,
    lab = false,
)

plot!(
    log.(ds),
    log.(1 ./ ds .^ 3),
    linewidth = 2,
    color = RGB(50 / 255, 50 / 255, 50 / 255),
    line = :dash,
    lab = false,
)

annotate!(
    x_min + (0.95) * (x_max - x_min),
    y_min + (0.95) * (y_max - y_min),
    text("(a)", 18),
)

savefig("plotters/plots/monoatomic/Monoatomic_T0_M-M.pdf")

# Finite T
ExD = readdlm("data/monoatomic_chain/ExD_M-M_Finite_T.dat")
QFT = readdlm("data/monoatomic_chain/QFT_M-M_Finite_T.dat")

plot(
    xaxis = (L"\ln d", font(14, "Serif")),
    yaxis = (L"\ln \left(E_I/E_I^1\right)", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    yticks = -10:5:0,
    xlims = (x_min, x_max),
    ylims = (y_min, y_max),
)
for ii = 1:length(Defect_Ms)
    plot!(
        log.(ds),
        log.(ExD[ii, :] ./ ExD[ii, 1]),
        color = colors[ii],
        label = false,
        linewidth = 2,
    )
    scatter!(
        log.(ds),
        log.(QFT[ii, :] ./ QFT[ii, 1]),
        color = colors[ii],
        label = latexstring("\$M = $(Defect_Ms[ii])\$"),
    )
end

plot!(
    log.(ds),
    log.(1 ./ ds),
    linewidth = 2,
    color = RGB(50 / 255, 50 / 255, 50 / 255),
    line = :dash,
    lab = false,
)

plot!(
    log.(ds),
    log.(1 ./ ds .^ 3),
    linewidth = 2,
    color = RGB(50 / 255, 50 / 255, 50 / 255),
    line = :dash,
    lab = false,
)

annotate!(
    x_min + (0.95) * (x_max - x_min),
    y_min + (0.95) * (y_max - y_min),
    text("(b)", 18),
)

savefig("plotters/plots/monoatomic/Monoatomic_Finite_T_M-M.pdf")

## Well-well interaction plots
x_min = -0.1;
x_max = log(maximum(ds)) + 0.1;
y_min = -5
y_max = 0.25

# T = 0
ExD = readdlm("data/monoatomic_chain/ExD_D-D_T0.dat")
QFT = readdlm("data/monoatomic_chain/QFT_D-D_T0.dat")


plot(
    xaxis = (L"\ln d", font(14, "Serif")),
    yaxis = (L"\ln \left(E_I/E_I^1\right)", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    yticks = -10:5:0,
    xlims = (x_min, x_max),
    ylims = (y_min, y_max),
)
for ii = 1:length(Defect_Δs)
    plot!(
        log.(ds),
        log.(ExD[ii, :] ./ ExD[ii, 1]),
        color = colors[ii],
        label = false,
        linewidth = 2,
    )
    scatter!(
        log.(ds),
        log.(QFT[ii, :] ./ QFT[ii, 1]),
        color = colors[ii],
        label = latexstring("\$\\Delta/k = $(Defect_Δs[ii])\$"),
    )
end

plot!(
    log.(ds),
    log.(1 ./ ds),
    linewidth = 2,
    color = RGB(50 / 255, 50 / 255, 50 / 255),
    line = :dash,
    lab = false,
)

annotate!(
    x_min + (0.95) * (x_max - x_min),
    y_min + (0.95) * (y_max - y_min),
    text("(a)", 18),
)

savefig("plotters/plots/monoatomic/Monoatomic_T0_D-D.pdf")

# Finite T
ExD = readdlm("data/monoatomic_chain/ExD_D-D_Finite_T.dat")
QFT = readdlm("data/monoatomic_chain/QFT_D-D_Finite_T.dat")

plot(
    xaxis = (L"\ln d", font(14, "Serif")),
    yaxis = (L"\ln \left(E_I/E_I^1\right)", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    yticks = -10:5:0,
    xlims = (x_min, x_max),
    ylims = (y_min, y_max),
)
for ii = 1:length(Defect_Ms)
    plot!(
        log.(ds),
        log.(ExD[ii, :] ./ ExD[ii, 1]),
        color = colors[ii],
        label = false,
        linewidth = 2,
    )
    scatter!(
        log.(ds),
        log.(QFT[ii, :] ./ QFT[ii, 1]),
        color = colors[ii],
        label = latexstring("\$\\Delta/k = $(Defect_Δs[ii])\$"),
    )
end

plot!(
    log.(ds),
    log.(1 ./ ds),
    linewidth = 2,
    color = RGB(50 / 255, 50 / 255, 50 / 255),
    line = :dash,
    lab = false,
)

annotate!(
    x_min + (0.95) * (x_max - x_min),
    y_min + (0.95) * (y_max - y_min),
    text("(b)", 18),
)

savefig("plotters/plots/monoatomic/Monoatomic_Finite_T_D-D.pdf")

## Mixed Impurities
x_min = -0.1;
x_max = log(maximum(ds)) + 0.1;
y_min = -10
y_max = 0.25

# T = 0
data = [
    readdlm("data/monoatomic_chain/D-D_T0_Mixed.dat"),
    readdlm("data/monoatomic_chain/M-D_T0_Mixed.dat"),
    readdlm("data/monoatomic_chain/M-M_T0_Mixed.dat"),
]

labels = ["\$ \\Delta\$-\$\\Delta\$", "\$M\$-\$\\Delta\$", "\$M\$-\$M\$"]

plot(
    xaxis = (L"\ln d", font(14, "Serif")),
    yaxis = (L"\ln \left(E_I/E_I^1\right)", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    yticks = -10:5:0,
    xlims = (x_min, x_max),
    ylims = (y_min, y_max),
)


for ii = 1:length(data)
    plot!(
        log.(ds),
        log.(data[ii][:, 1] ./ data[ii][1, 1]),
        color = colors[ii],
        label = false,
        linewidth = 2,
    )
    scatter!(
        log.(ds),
        log.(data[ii][:, 2] ./ data[ii][1, 2]),
        color = colors[ii],
        label = latexstring(labels[ii]),
    )
end

plot!(
    log.(ds),
    log.(1 ./ ds),
    linewidth = 2,
    color = RGB(50 / 255, 50 / 255, 50 / 255),
    line = :dash,
    lab = false,
)

plot!(
    log.(ds),
    log.(1 ./ ds .^ 3),
    linewidth = 2,
    color = RGB(50 / 255, 50 / 255, 50 / 255),
    line = :dash,
    lab = false,
)

annotate!(
    x_min + (0.95) * (x_max - x_min),
    y_min + (0.95) * (y_max - y_min),
    text("(a)", 18),
)

savefig("plotters/plots/monoatomic/Monoatomic_T0_Mixed.pdf")


# Finite T
data = [
    readdlm("data/monoatomic_chain/D-D_Finite_T_Mixed.dat"),
    readdlm("data/monoatomic_chain/M-D_Finite_T_Mixed.dat"),
    readdlm("data/monoatomic_chain/M-M_Finite_T_Mixed.dat"),
]

labels = ["\$ \\Delta\$-\$\\Delta\$", "\$M\$-\$\\Delta\$", "\$M\$-\$M\$"]

plot(
    xaxis = (L"\ln d", font(14, "Serif")),
    yaxis = (L"\ln \left(E_I/E_I^1\right)", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    yticks = -10:5:0,
    xlims = (x_min, x_max),
    ylims = (y_min, y_max),
)


for ii = 1:length(data)
    plot!(
        log.(ds),
        log.(data[ii][:, 1] ./ data[ii][1, 1]),
        color = colors[ii],
        label = false,
        linewidth = 2,
    )
    scatter!(
        log.(ds),
        log.(data[ii][:, 2] ./ data[ii][1, 2]),
        color = colors[ii],
        label = latexstring(labels[ii]),
    )
end

plot!(
    log.(ds),
    log.(1 ./ ds),
    linewidth = 2,
    color = RGB(50 / 255, 50 / 255, 50 / 255),
    line = :dash,
    lab = false,
)

plot!(
    log.(ds),
    log.(1 ./ ds .^ 3),
    linewidth = 2,
    color = RGB(50 / 255, 50 / 255, 50 / 255),
    line = :dash,
    lab = false,
)

annotate!(
    x_min + (0.95) * (x_max - x_min),
    y_min + (0.95) * (y_max - y_min),
    text("(b)", 18),
)

savefig("plotters/plots/monoatomic/Monoatomic_Finite_T_Mixed.pdf")

## Cluster Plot
x_min = -0.1;
x_max = maximum(Ds) - minimum(Ds) + 0.1;
y_min = -7.5
y_max = 4
data = readdlm("data/monoatomic_chain/Cluster.dat")
plot(
    xaxis = (L"d", font(14, "Serif")),
    yaxis = (L"\left(E_I - E_I^0\right)\times 10^4", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    xlims = (x_min, x_max),
    ylims = (y_min, y_max),
)
for ii = 1:length(Ms_cluster)
    plot!(
        Ds .- minimum(Ds),
        (data[:, ii] .- data[1, ii]) * 1e4,
        color = colors[ii],
        linewidth = 2,
        label = L"M = %$(Ms_cluster[ii])",
    )

end
annotate!(
    x_min + (0.05) * (x_max - x_min),
    y_min + (0.95) * (y_max - y_min),
    text("(c)", 18),
)

savefig("plotters/plots/monoatomic/Monoatomic_cluster.pdf")

## Heatmap
ΔM_heatmap_data = readdlm("data/monoatomic_chain/M-Delta-Heatmap.dat")
ΔΔ_heatmap_data = readdlm("data/monoatomic_chain/Delta-Delta-Heatmap.dat")

heatmap(
    MΔ,
    MΔ,
    ΔΔ_heatmap_data .* 1e2,
    color = reverse(cgrad(:Blues_9)),
    xaxis = (L"\Delta_L", font(14, "Serif")),
    yaxis = (L"\Delta_R", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    colorbar_title = (L"$(E_I/\Omega)\times 10^2$"),
    colorbar_titlefontsize = 18
)

annotate!(
    (0.95) * maximum(MΔ),
    (0.95) * maximum(MΔ),
    text("(a)", 18, color = :white),
)

savefig("plotters/plots/monoatomic/D-D_heatmap.pdf")

heatmap(
    MΔ,
    MΔ,
    ΔM_heatmap_data .* 1e2,
    color = :RdBu_7,
    xaxis = (L"M", font(14, "Serif")),
    yaxis = (L"\Delta", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    colorbar_title = (L"$(E_I/\Omega)\times 10^2$"),

)
annotate!(
    (0.95) * maximum(MΔ),
    (0.95) * maximum(MΔ),
    text("(b)", 18, color = :white),
)

savefig("plotters/plots/monoatomic/D-M_heatmap.pdf")
