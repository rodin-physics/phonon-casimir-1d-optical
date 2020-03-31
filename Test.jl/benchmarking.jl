include("../src/general.jl")

μ = 1.0
M = 2.0
N = 500

T = 1e-12

s = System([μ, 2 * μ], [Impurity(1, 1, M), Impurity(10, 1, 2 * M)], T, N)

@code_warntype exact_F(s)
@time exact_F(s)


@time Ξ(0.1 + 1im * η, s, 1, 1)
@code_warntype Ξ(0.1 + 1im * η, s, 1, 1)

#
# θs = range(-π, stop = π; length = 100)
#
# aco = map(x -> modes([1, 1.5], x).values[1], θs)
# opt = map(x -> modes([1, 1.5], x).values[2], θs)
# pyplot()
#
# plot(θs, sqrt.(aco))
# plot!(θs, sqrt.(opt))
# savefig("Test.pdf")
#
# @code_warntype modes([2, 1], 0.4)
# F = modes([2, 1], 0.4)
#
# @time Q2_sθ([2, 1], 1, 0.4)
#
# Imps = [Impurity(1, 1, M), Impurity(3, 1, M)]
# M_List = repeat([1], 10)
#
#
# M_List = reduce(x -> M_List[nAtom*(x.pos-1)+x.n] = x.λ, Imps)
#
#
#
# # Replace the pristine masses by the impurities
# for ii in Imps
#     coord = 1 * (ii.pos - 1) + ii.n
#     M_List[coord] = ii.λ
# end
# M_List
# # M_List = reduce(x -> M_List[nAtom*(x.pos-1)+x.n] = x.λ, Imps; init = M_List)
