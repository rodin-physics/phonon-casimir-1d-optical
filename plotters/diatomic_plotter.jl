include("../calculations/diatomic_calc_settings.jl")
pgfplotsx()

## Cluster Plot
x_min = -0.1;
x_max = maximum(Ds) - minimum(Ds) + 0.1;
y_min = -1
y_max = 0.5
data = readdlm("data/diatomic_chain/Cluster.dat")
plot(
    xaxis = (L"d", font(14, "Serif")),
    yaxis = (L"\left(E_I - E_I^0\right)/\Omega\times 10^3", font(14, "Serif")),
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
        (data[:, ii] .- data[1, ii]) * 1e3,
        color = colors[ii],
        linewidth = 2,
        label = L"M = %$(Ms_cluster[ii])",
    )

end
annotate!(
    x_min + (0.95) * (x_max - x_min),
    y_min + (0.95) * (y_max - y_min),
    text("(c)", 18),
)
savefig("plotters/plots/diatomic/Diatomic_cluster.pdf")

## Mass-mass interaction plots

heavy_light_T0 = readdlm("data/diatomic_chain/Heavy-Light_T0.dat")
light_heavy_T0 = readdlm("data/diatomic_chain/Light-Heavy_T0.dat")

heavy_light_Finite_T = readdlm("data/diatomic_chain/Heavy-Light_Finite_T.dat")
light_heavy_Finite_T = readdlm("data/diatomic_chain/Light-Heavy_Finite_T.dat")

x_min = -0.1;
x_max = log(maximum(ds)) + 0.1;
y_min = -.25
y_max = 1.1

plot(
    xaxis = (L"\ln d", font(14, "Serif")),
    yaxis = (L"E_I/E_I^1", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    yticks = 0:1,
    legend = :topright,
    xlims = (x_min, x_max),
    ylims = (y_min, y_max),
)

plot!(
    log.(ds),
    heavy_light_T0[:, 1] ./ heavy_light_T0[1, 1],
    linewidth = 2,
    color = colors[5],
    lab = false,
)
plot!(
    log.(ds),
    light_heavy_T0[:, 1] ./ light_heavy_T0[1, 1],
    linewidth = 2,
    color = colors[4],
    lab = false,
)
scatter!(
    log.(ds),
    heavy_light_T0[:, 2] ./ heavy_light_T0[1, 2],
    color = colors[5],
    lab = "fix heavy",
)
scatter!(
    log.(ds),
    light_heavy_T0[:, 2] ./ light_heavy_T0[1, 2],
    color = colors[4],
    lab = "fix light",
)

savefig("plotters/plots/diatomic/Diatomic_T0_M-M_Linear.pdf")

## Monotonic interaction plots
x_min = -0.1;
x_max = log(maximum(ds / 2)) + 0.1;
y_min = -12
y_max = 0.1
# T = 0
plot(
    xaxis = (L"\ln d_u", font(14, "Serif")),
    yaxis = (L"\ln \left(E_I/E_I^1\right)", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    yticks = -10:5:0,
    xlims = (x_min, x_max),
    ylims = (y_min, y_max),
)
plot!(
    log.(1:length(ds)/2),
    log.(heavy_light_T0[1:2:end, 1] ./ heavy_light_T0[1, 1]),
    linewidth = 2,
    color = colors[1],
    lab = false,
)

plot!(
    log.(1:length(ds)/2),
    log.(light_heavy_T0[2:2:end, 1] ./ light_heavy_T0[2, 1]),
    linewidth = 2,
    color = colors[2],
    lab = false,
)
plot!(
    log.(1:length(ds)/2),
    log.(heavy_light_T0[2:2:end, 1] ./ heavy_light_T0[2, 1]),
    linewidth = 2,
    color = colors[3],
    lab = false,
)

scatter!(
    log.(1:length(ds)/2),
    log.(heavy_light_T0[1:2:end, 1] ./ heavy_light_T0[1, 1]),
    color = colors[1],
    lab = "heavy-heavy",
)

scatter!(
    log.(1:length(ds)/2),
    log.(light_heavy_T0[2:2:end, 1] ./ light_heavy_T0[2, 1]),
    color = colors[2],
    lab = "light-light",
)
scatter!(
    log.(1:length(ds)/2),
    log.(heavy_light_T0[2:2:end, 1] ./ heavy_light_T0[2, 1]),
    color = colors[3],
    lab = "heavy-light",
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

savefig("plotters/plots/diatomic/Diatomic_Monotonic_Interaction_T0.pdf")

# Finite T
plot(
    xaxis = (L"\ln d_u", font(14, "Serif")),
    yaxis = (L"\ln \left(E_I/E_I^1\right)", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    yticks = -10:5:0,
    xlims = (x_min, x_max),
    ylims = (y_min, y_max),
)
plot!(
    log.(1:length(ds)/2),
    log.(heavy_light_Finite_T[1:2:end, 1] ./ heavy_light_Finite_T[1, 1]),
    linewidth = 2,
    color = colors[1],
    lab = false,
)

plot!(
    log.(1:length(ds)/2),
    log.(light_heavy_Finite_T[2:2:end, 1] ./ light_heavy_Finite_T[2, 1]),
    linewidth = 2,
    color = colors[2],
    lab = false,
)
plot!(
    log.(1:length(ds)/2),
    log.(heavy_light_Finite_T[2:2:end, 1] ./ heavy_light_Finite_T[2, 1]),
    linewidth = 2,
    color = colors[3],
    lab = false,
)

scatter!(
    log.(1:length(ds)/2),
    log.(heavy_light_Finite_T[1:2:end, 1] ./ heavy_light_Finite_T[1, 1]),
    color = colors[1],
    lab = "heavy-heavy",
)

scatter!(
    log.(1:length(ds)/2),
    log.(light_heavy_Finite_T[2:2:end, 1] ./ light_heavy_Finite_T[2, 1]),
    color = colors[2],
    lab = "light-light",
)
scatter!(
    log.(1:length(ds)/2),
    log.(heavy_light_Finite_T[2:2:end, 1] ./ heavy_light_Finite_T[2, 1]),
    color = colors[3],
    lab = "heavy-light",
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

savefig("plotters/plots/diatomic/Diatomic_Monotonic_Interaction_Finite_T.pdf")
