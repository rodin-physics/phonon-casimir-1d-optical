using LinearAlgebra
using Plots

function pristineRowDisp(Ms, N, Ωx)
    dv = 2 .* ones(N)
    ev = -ones(N - 1)
    U_Mat = SymTridiagonal(dv, ev) + zeros(N, N)
    # make system periodic
    U_Mat[1, N] = -1
    U_Mat[N, 1] = -1
    # prepare the matrix of masses
    M_Mat = repeat([Ms], N) |> Diagonal

    # Calculate the eigenvalues of the system which give
    # the squares of the energies
    Ω2 = eigvals(inv(M_Mat) * U_Mat)
    Ω = sqrt.(abs.(Ω2)) # in units √(k / μ)
    # plotly plots in browser
    # plotly()
    # plot(res[1:N])
    # plot!(reverse(res[N+1:2*N]))
    plot(1:N, Ω)
end

pristineRowDisp(1, 3, 2)
