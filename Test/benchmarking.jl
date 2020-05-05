include("../src/general.jl")

μ = 1.0
M = 3.0
N = 500

T = 1e-12
s = System([μ, μ], [Impurity(1, 1, M), Impurity(3, 2, M)], T, N)
sFar = System([μ, μ], [Impurity(1, 1, M), Impurity(250, 1, M)], T, N)

# About 3e-5 seconds
@time modes(s.Ms, 1/1.2)

# Same as above
z = 1.2 + 1im * η
@time Ξ_jl_integrand(z, s.Ms, 1.22, 1, 1)

# About 1e-2 seconds
@time Ξ_jl(z, s, 1, 2)

# Should be about 4 times longer than above
@time Ξ(z, s)

# Should be close to above
@time F_I_Integrand_T0(z, s)

# Sub-second for short distances. Order of seconds for larger
@time F_I_T0(s)

# Depends on the size of the system. Order of seconds
@time exact_F(s) - exact_F(sFar)
