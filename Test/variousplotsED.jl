# include("monoatomic_library.jl")
include("../src/exactDiag1D.jl")
include("../src/general.jl")

# Comparing for the 0 temperature case

Ds = 2:50
Ds_Exact = 2 : 3 : 50   # Distances used for the exact diagonalization
N = 1000                # Chain length for the exact diagonalization
# T in units of Ω0
T = 1e-12;
Ds = Ds_Exact
# Impurity masses in units of m (chain masses)
Ms = [1/3, 4/3, 4, 10, 1000]
diMs = [1., 3., 12., 50.]


# Bounding power laws
r_1 = 1 ./ Ds
r_3 = 1 ./ (Ds.^3)

# Comparing for the 0 temperature case

Ds = 2:50
Ds_Exact = 2 : 3 : 50   # Distances used for the exact diagonalization
N = 1000                # Chain length for the exact diagonalization
# T in units of Ω0
T = 1e-12;
Ds = Ds_Exact
# Impurity masses in units of m (chain masses)
Ms = [1/3, 4/3, 4, 10, 1000]

# Bounding power laws
r_1 = 1 ./ Ds
r_3 = 1 ./ (Ds.^3)

## Plotting
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

function two_diImps(D, impM)
    Imps = [Impurity(1, 1, impM), Impurity(1 + floor(Int, D/2), 1 + Int(D % 2), impM)]
    return Imps
end

jjkkii = Array{Array{Float64}}(undef, 4, 4, 5)
fill!(jjkkii, [])


@time begin
for jj = 1:length(diMs)
    println(jj)
    for kk = jj:length(diMs)
        for ii = 1:(length(Ms)- 0)
            println(ii)
            diM = [diMs[jj], diMs[kk]];
            M = Ms[ii];
            r =  map(x -> exact_F(System(diM, two_diImps(x, M), T, N)), Ds_Exact);
            jjkkii[jj, kk, ii] = r
            # Plots.display(plot!(log.(Ds_Exact), real(log.(complex(r))),
            #     color = colors[ii],
            #     lab = "",
            #     markershape = :cross
            #     ))
        end
    end
end
println("done!")
end


plot(
    xaxis = (L"$\ln\, D$", font(14, "Serif")),
    yaxis = (L"$\ln \left(F_I/F_I^1\right) $", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    yticks = -15 : 5 : 0,
    ylims = (-20, 5)
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

plotData = jjkkii[1, 1, :]
for ii = 1:(length(Ms) - 0)
    Plots.display(plot!(log.(Ds_Exact), real(log.(complex(plotData[ii]))),
        color = colors[ii],
        lab = "",
        markershape = :cross
        ))
end

# Plots.display(plot!(log.(Ds_Exact), testArr,
#     color = colors[1],
#     lab = "",
#     markershape = :cross
#     ))
for jj = 1:(length(diMs) - 0)
    println(jj)
    for kk = jj:(length(diMs) - 0)
        plot(
            xaxis = (L"$\ln\, D$", font(14, "Serif")),
            yaxis = (L"$\ln \left(F_I/F_I^1\right) $", font(14, "Serif")),
            xtickfont = font(12, "Serif"),
            ytickfont = font(12, "Serif"),
            legendfont = font(12, "Serif"),
            legend = :bottomleft,
            yticks = -15 : 5 : 0,
            ylims = (-20, 5)
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

        plotData = jjkkii[jj, kk, :]
        for ii = 1:(length(Ms) - 0)
            Plots.display(plot!(log.(Ds_Exact), real(log.(complex(plotData[ii]))),
                color = colors[ii],
                lab = "",
                markershape = :cross
                ))
        end
        savefig("plotEDv2" * "_" * string(jj) * "_" * string(kk) * "Ms.png")
    end
end
println("done!")

plotData = jjkkii[1, 1, :]
for ii in 1:5
    plot(log.(Ds_Exact), real(log.(complex(plotData[ii]))))
end
plot(log.(Ds_Exact), real(log.(complex(plotData[5]))))
savefig("plot3_3_Ms.png")

testArr = real(log.(complex(jjkkii[1, 4, 3])))

b = 2
if b == 2
    println("hi")
end
