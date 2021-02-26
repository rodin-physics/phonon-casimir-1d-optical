using Distributed
## Set the number of Procs
nProcs = 8;
if nprocs() < nProcs
    addprocs(nProcs - nprocs())
end
@everywhere include("calculations/monoatomic_calc_settings.jl")

## Pair interaction
# Mass-Mass
# T = 0
res_ExD = @showprogress pmap(Defect_Ms_mat, ds_M_mat) do M, d
    Imp_1 = Impurity(1, 1, M, 0)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), M, 0)
    return exact_F_I(System(Ms, [Imp_1, Imp_2], 0, N, K))
end
writedlm("data/monoatomic_chain/ExD_M-M_T0.dat", res_ExD)

res_QFT = @showprogress pmap(Defect_Ms_mat, ds_M_mat) do M, d
    Imp_1 = Impurity(1, 1, M, 0)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), M, 0)
    return E_I(System(Ms, [Imp_1, Imp_2], 0, N, K))
end
writedlm("data/monoatomic_chain/QFT_M-M_T0.dat", res_QFT)

# Finite T
res_ExD = @showprogress pmap(Defect_Ms_mat, ds_M_mat) do M, d
    Imp_1 = Impurity(1, 1, M, 0)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), M, 0)
    return exact_F_I(System(Ms, [Imp_1, Imp_2], T, N, K))
end
writedlm("data/monoatomic_chain/ExD_M-M_Finite_T.dat", res_ExD)

res_QFT = @showprogress pmap(Defect_Ms_mat, ds_M_mat) do M, d
    Imp_1 = Impurity(1, 1, M, 0)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), M, 0)
    return E_I(System(Ms, [Imp_1, Imp_2], T, N, K))
end
writedlm("data/monoatomic_chain/QFT_M-M_Finite_T.dat", res_QFT)

# Well-Well
# T = 0
res_ExD = @showprogress pmap(Defect_Δs_mat, ds_Δ_mat) do Δ, d
    Imp_1 = Impurity(1, 1, m, Δ)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), m, Δ)
    return (
        exact_F_I(System(Ms, [Imp_1, Imp_2], 0, N, K)) -
        exact_neg_TS_I(System(Ms, [Imp_1, Imp_2], 0, N, K))
    )
end
writedlm("data/monoatomic_chain/ExD_D-D_T0.dat", res_ExD)

res_QFT = @showprogress pmap(Defect_Δs_mat, ds_M_mat) do Δ, d
    Imp_1 = Impurity(1, 1, m, Δ)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), m, Δ)
    return E_I(System(Ms, [Imp_1, Imp_2], 0, N, K))
end
writedlm("data/monoatomic_chain/QFT_D-D_T0.dat", res_QFT)

# Finite T
res_ExD = @showprogress pmap(Defect_Δs_mat, ds_Δ_mat) do Δ, d
    Imp_1 = Impurity(1, 1, m, Δ)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), m, Δ)
    return (
        exact_F_I(System(Ms, [Imp_1, Imp_2], T, N, K)) -
        exact_neg_TS_I(System(Ms, [Imp_1, Imp_2], T, N, K))
    )
end
writedlm("data/monoatomic_chain/ExD_D-D_Finite_T.dat", res_ExD)

res_QFT = @showprogress pmap(Defect_Δs_mat, ds_M_mat) do Δ, d
    Imp_1 = Impurity(1, 1, m, Δ)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), m, Δ)
    return E_I(System(Ms, [Imp_1, Imp_2], T, N, K))
end
writedlm("data/monoatomic_chain/QFT_D-D_Finite_T.dat", res_QFT)

## Mixed configurations
# T = 0
res_Mass_Mass = @showprogress pmap(ds) do d
    Imp_1 = Impurity(1, 1, Imp_M, 0)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), Imp_M, 0)
    return [
        exact_F_I(System(Ms, [Imp_1, Imp_2], 0, N, K)) -
        exact_neg_TS_I(System(Ms, [Imp_1, Imp_2], 0, N, K)),
        E_I(System(Ms, [Imp_1, Imp_2], 0, N, K)),
    ]
end
writedlm("data/monoatomic_chain/M-M_T0_Mixed.dat", res_Mass_Mass)

res_Mass_Delta = @showprogress pmap(ds) do d
    Imp_1 = Impurity(1, 1, Imp_M, 0)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), m, Imp_Δ)
    return [
        exact_F_I(System(Ms, [Imp_1, Imp_2], 0, N, K)) -
        exact_neg_TS_I(System(Ms, [Imp_1, Imp_2], 0, N, K)),
        E_I(System(Ms, [Imp_1, Imp_2], 0, N, K)),
    ]
end
writedlm("data/monoatomic_chain/M-D_T0_Mixed.dat", res_Mass_Delta)

res_Delta_Delta = @showprogress pmap(ds) do d
    Imp_1 = Impurity(1, 1, m, Imp_Δ)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), m, Imp_Δ)
    return [
        exact_F_I(System(Ms, [Imp_1, Imp_2], 0, N, K)) -
        exact_neg_TS_I(System(Ms, [Imp_1, Imp_2], 0, N, K)),
        E_I(System(Ms, [Imp_1, Imp_2], 0, N, K)),
    ]
end
writedlm("data/monoatomic_chain/D-D_T0_Mixed.dat", res_Delta_Delta)

# Finite T
res_Mass_Mass = @showprogress pmap(ds) do d
    Imp_1 = Impurity(1, 1, Imp_M, 0)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), Imp_M, 0)
    return [
        exact_F_I(System(Ms, [Imp_1, Imp_2], T, N, K)) -
        exact_neg_TS_I(System(Ms, [Imp_1, Imp_2], T, N, K)),
        E_I(System(Ms, [Imp_1, Imp_2], T, N, K)),
    ]
end
writedlm("data/monoatomic_chain/M-M_Finite_T_Mixed.dat", res_Mass_Mass)


res_Mass_Delta = @showprogress pmap(ds) do d
    Imp_1 = Impurity(1, 1, Imp_M, 0)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), m, Imp_Δ)
    return [
        exact_F_I(System(Ms, [Imp_1, Imp_2], T, N, K)) -
        exact_neg_TS_I(System(Ms, [Imp_1, Imp_2], T, N, K)),
        E_I(System(Ms, [Imp_1, Imp_2], T, N, K)),
    ]
end
writedlm("data/monoatomic_chain/M-D_Finite_T_Mixed.dat", res_Mass_Delta)

res_Delta_Delta = @showprogress pmap(ds) do d
    Imp_1 = Impurity(1, 1, m, Imp_Δ)
    Imp_2 = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), m, Imp_Δ)
    return [
        exact_F_I(System(Ms, [Imp_1, Imp_2], T, N, K)) -
        exact_neg_TS_I(System(Ms, [Imp_1, Imp_2], T, N, K)),
        E_I(System(Ms, [Imp_1, Imp_2], T, N, K)),
    ]
end
writedlm("data/monoatomic_chain/D-D_Finite_T_Mixed.dat", res_Delta_Delta)

## Cluster
res = @showprogress pmap(Ds_cluster_mat, Ms_cluster_mat) do d, m
    middle_imp = Impurity(1 + Integer(floor(d / 2)), 1 + isodd(d), m, 0)
    E_I(System(Ms, [left_imp, right_imp, middle_imp], 0, N, K))
end
writedlm("data/monoatomic_chain/Cluster.dat", res)

## Heatmap

res_ΔΔ = @showprogress pmap(MΔ_mat, MΔ_mat') do Δ1, Δ2
    Imp_1 = Impurity(1, 1, m, Δ1)
    Imp_2 = Impurity(1, 2, m, Δ2)
    return E_I(System(Ms, [Imp_1, Imp_2], 0, N, K))

end

writedlm("data/monoatomic_chain/Delta-Delta-Heatmap.dat", res_ΔΔ)


res_ΔM = @showprogress pmap(MΔ_mat, MΔ_mat') do Δ, M
    Imp_1 = Impurity(1, 1, m, Δ)
    Imp_2 = Impurity(1, 2, M, 0)
    return E_I(System(Ms, [Imp_1, Imp_2], 0, N, K))

end

writedlm("data/monoatomic_chain/M-Delta-Heatmap.dat", res_ΔM)
