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
jj = 4
kk = 4
imax = length(Ms)
# if kk >= 3
#     imax = length(Ms) - 2
# end
## kk is the heavier atom
for ii = 1:(imax - 0)
    println(ii)
    diM = [diMs[jj], diMs[kk]];
    M = Ms[ii];
    # F_halfway = F_I_T0(System(diM, two_diImps(floor(Int, N/2), M), T, N))
    F_halfway = 0
    F0 = F_I_T0(System(diM, two_diImps(1, M), T, N)) - F_halfway
    r =  map(x -> F_I_T0(System(diM, two_diImps(x, M), T, N)) - F_halfway, Ds_Exact) ./ F0;
    jjkkii[jj, kk, ii] = r
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

plotData = jjkkii[4, 4, :]
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

savefig("plot50_50_Ms.png")

testArr = real(log.(complex(jjkkii[1, 4, 3])))

b = 2
if b == 2
    println("hi")
end
