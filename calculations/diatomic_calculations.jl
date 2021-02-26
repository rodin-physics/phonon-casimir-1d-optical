using Distributed
## Set the number of Procs
nProcs = 8;
if nprocs() < nProcs
    addprocs(nProcs - nprocs())
end
@everywhere include("calculations/diatomic_calc_settings.jl")

## Pair interaction
# T = 0
res_light_heavy = @showprogress pmap(ds) do d
    Imp_1 = Impurity(1, 1, Imp_M, 0)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), Imp_M, 0)
    return [
        exact_F_I(System([m1, m2], [Imp_1, Imp_2], 0, N, K)),
        E_I(System([m1, m2], [Imp_1, Imp_2], 0, N, K)),
    ]
end

res_heavy_light = @showprogress pmap(ds) do d
    Imp_1 = Impurity(1, 1, Imp_M, 0)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), Imp_M, 0)
    return [
        exact_F_I(System([m2, m1], [Imp_1, Imp_2], 0, N, K)),
        E_I(System([m2, m1], [Imp_1, Imp_2], 0, N, K)),
    ]
end
writedlm("data/diatomic_chain/Light-Heavy_T0.dat", res_light_heavy)
writedlm("data/diatomic_chain/Heavy-Light_T0.dat", res_heavy_light)

# Finite T
res_light_heavy = @showprogress pmap(ds) do d
    Imp_1 = Impurity(1, 1, Imp_M, 0)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), Imp_M, 0)
    return [
        exact_F_I(System([m1, m2], [Imp_1, Imp_2], T, N, K)),
        E_I(System([m1, m2], [Imp_1, Imp_2], T, N, K)),
    ]
end

res_heavy_light = @showprogress pmap(ds) do d
    Imp_1 = Impurity(1, 1, Imp_M, 0)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), Imp_M, 0)
    return [
        exact_F_I(System([m2, m1], [Imp_1, Imp_2], T, N, K)),
        E_I(System([m2, m1], [Imp_1, Imp_2], T, N, K)),
    ]
end
writedlm("data/diatomic_chain/Light-Heavy_Finite_T.dat", res_light_heavy)
writedlm("data/diatomic_chain/Heavy-Light_Finite_T.dat", res_heavy_light)



## Cluster
res = @showprogress pmap(Ds_cluster_mat, Ms_cluster_mat) do d, m
    middle_imp = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), m, 0)
    E_I(System(Ms, [left_imp, right_imp, middle_imp], 0, N, K))
end
writedlm("data/diatomic_chain/Cluster.dat", res)
