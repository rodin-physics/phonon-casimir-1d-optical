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
diMs = [1, 3, 12, 50]

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
jjkkiiED = Array{Array{Float64}}(undef, 4, 4, 5)
fill!(jjkkiiED, [])

@time begin
for jj = 2:3
    print("jj", jj)
    for kk = jj:4
        print("kk", kk)
        imax = length(Ms)
        ## kk is the heavier atom
        for ii = 1:(imax - 0)
            println("ii", ii)
            diM = [diMs[jj], diMs[kk]];
            M = Ms[ii];

            # QFT
            F0 = F_I_T0(System(diM, two_diImps(1, M), T, N))
            r =  map(x -> F_I_T0(System(diM, two_diImps(x, M), T, N)), Ds_Exact) ./ F0;
            jjkkii[jj, kk, ii] = r
            println("QFT Done")
            # ED exact_FI(s
            # F0ED = exact_FI(System(diM, two_diImps(1, M), T, N))
            # rED =  map(x -> exact_FI(System(diM, two_diImps(x, M), T, N)), Ds_Exact) ./ F0ED;
            # jjkkiiED[jj, kk, ii] = rED
            # println("ED Done")
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

@time begin
jj = 2
kk = 3
imax = length(Ms)
## kk is the heavier atom
for ii = 3:(imax - 0)
    println("ii", ii)
    diM = [diMs[jj], diMs[kk]];
    M = Ms[ii];

    # QFT
    F0 = F_I_T0(System(diM, two_diImps(1, M), T, N))
    r =  map(x -> F_I_T0(System(diM, two_diImps(x, M), T, N)), Ds_Exact) ./ F0;
    jjkkii[jj, kk, ii] = r
    println("QFT Done")
    # ED exact_FI(s
    # F0ED = exact_FI(System(diM, two_diImps(1, M), T, N))
    # rED =  map(x -> exact_FI(System(diM, two_diImps(x, M), T, N)), Ds_Exact) ./ F0ED;
    # jjkkiiED[jj, kk, ii] = rED
    # println("ED Done")
    # Plots.display(plot!(log.(Ds_Exact), real(log.(complex(r))),
    #     color = colors[ii],
    #     lab = "",
    #     markershape = :cross
    #     ))
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

jj = 2
kk = 4
plotDataQFT = jjkkii[jj, kk, :]
plotDataED = jjkkiiED[jj, kk, :]
for ii = 2:(length(Ms) - 1)
    Plots.display(plot!(log.(Ds_Exact), real(log.(complex(plotDataQFT[ii]))),
        color = colors[ii],
        lab = "",
        markershape = :cross
        ))
end

for ii = 2:(length(Ms) - 1)
    Plots.display(plot!(log.(Ds_Exact), real(log.(complex(plotDataED[ii]))),
        color = colors[ii],
        lab = "",
        markershape = :hexagon
        ))
end

# Plots.display(plot!(log.(Ds_Exact), testArr,
#     color = colors[1],
#     lab = "",
#     markershape = :cross
#     ))


# Ms = [1/3, 4/3, 4, 10, 1000]
# diMs = [1, 3, 12, 50]
savefig("plot_diatomic_EDvQFT_4o3_10.png")


using DelimitedFiles

writedlm( "diatomic_QFT_4o3_10_Ds=2-3-50.csv",  plotDataQFT, ',')
writedlm( "diatomic_ED_4o3_10_Ds=2-3-50.csv",  plotDataED, ',')
writedlm( "diatomic_jjkkii_tentative",  jjkkii, ',')


testArr = real(log.(complex(jjkkii[1, 4, 3])))

b = 2
if b == 2
    println("hi")
end
