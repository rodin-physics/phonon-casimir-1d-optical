using LinearAlgebra
using Plots

function dispersion_plot(Ms, N)
    # Prepare a pristine chain potential energy matrix
    dv = 2 .* ones(2N)
    ev = -ones(2N - 1)
    U_Mat = SymTridiagonal(dv, ev) + zeros(2N, 2N)
    # Make the system periodic by coupling the first and the last mass
    U_Mat[1, 2N] = -1
    U_Mat[2N, 1] = -1
    # Prepare the matrix of masses
    M_Mat = repeat(Ms, N) |> Diagonal

    # Calculate the eigenvalues of the system which give
    # the squares of the energies
    Ω2 = eigvals(inv(M_Mat) * U_Mat)
    Ω = sqrt.(abs.(Ω2)) # in units √(k / μ)
    # plotly plots in browser
    # plotly()
    # plot(res[1:N])
    # plot!(reverse(res[N+1:2*N]))
    plot(Ω[1:N])
    plot!(reverse(Ω[N+1:2*N]))
end

function exact_F(Ms, N, Imps, T)
    # Prepare a pristine chain potential energy matrix
    dv = 2 .* ones(2N)
    ev = -ones(2N - 1)
    U_Mat = SymTridiagonal(dv, ev) + zeros(2N, 2N)
    # Make the system periodic by coupling the first and the last mass
    U_Mat[1, 2N] = -1
    U_Mat[2N, 1] = -1
    # Prepare the matrix of masses
    M_Mat = repeat(Ms, N)
    # Replace the pristine masses by the impurities
    for ii in Imps
        coord = length(Ms) * (ii.pos - 1) + ii.n
        M_Mat[coord] = ii.λ
    end
    M_Mat = Diagonal(M_Mat)
    # Calculate the eigenvalues of the system which give
    # the squares of the energies
    Ω2 = eigvals(inv(M_Mat) * U_Mat)
    Ω = sqrt.(abs.(Ω2)) # in units √(k / μ)
    # The free energy for each mode consists of the vacuum portion Ω / 2 and
    # the finite-T portion T * log(1 - exp(-Ω / T))
    filter!(x -> x > 1e-12, Ω)

    total_energy = T * sum(log.(1 .- exp.(-Ω ./ T))) + sum(Ω) / 2
    return total_energy
end

# nPts = 400 # Number of unit cells
# Ms = [1, 1]
#
# dispersion_plot(Ms, nPts)
