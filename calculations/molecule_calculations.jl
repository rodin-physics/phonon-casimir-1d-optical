using Distributed

## Set the number of Procs
nProcs = 4;
if nprocs() < nProcs
    addprocs(nProcs - nprocs())
end

@everywhere include("calculations/molecule_calc_settings.jl")
## Calculations

# Zero temperature
# ExD
res_F_I_exact_T0_M = @showprogress pmap(
    (K, M) -> F_I_exact_molecule(0, M, K, 0),
    Ks_exact_mat,
    Ms_exact_mat,
)
writedlm("data/molecule/F_I_exact_T0_M.dat", res_F_I_exact_T0_M)

res_F_I_exact_T0_log_M = @showprogress pmap(
    (K, M) -> F_I_exact_molecule(0, M, K, 0),
    Ks_exact_mat,
    log_Ms_exact_mat,
)
writedlm("data/molecule/F_I_exact_T0_log_M.dat", res_F_I_exact_T0_log_M)

res_F_I_exact_T0_δ = @showprogress pmap(
    (K, δ) -> F_I_exact_molecule(δ, 1, K, 0),
    Ks_exact_mat,
    δs_exact_mat,
)
writedlm("data/molecule/F_I_exact_T0_delta.dat", res_F_I_exact_T0_δ)

#QFT
res_F_I_QFT_T0_M = @showprogress pmap(
    (K, M) -> F_I_molecule(0, M, K, 0),
    Ks_QFT_mat,
    Ms_QFT_mat,
)
writedlm("data/molecule/F_I_QFT_T0_M.dat", res_F_I_QFT_T0_M)

res_F_I_QFT_T0_log_M = @showprogress pmap(
    (K, M) -> F_I_molecule(0, M, K, 0),
    Ks_QFT_mat,
    log_Ms_QFT_mat,
)
writedlm("data/molecule/F_I_QFT_T0_log_M.dat", res_F_I_QFT_T0_log_M)

res_F_I_QFT_T0_δ = @showprogress pmap(
    (K, δ) -> F_I_molecule(δ, 1, K, 0),
    Ks_QFT_mat,
    δs_QFT_mat,
)
writedlm("data/molecule/F_I_QFT_T0_delta.dat", res_F_I_QFT_T0_δ)

# Function of temperature
# ExD
res_F_I_exact_T_single_M = @showprogress pmap(
    (K, T) -> F_I_exact_molecule(0, single_M, K, T),
    Ks_exact_mat,
    Ts_exact_mat,
)
writedlm("data/molecule/F_I_exact_T_single_M.dat", res_F_I_exact_T_single_M)

res_F_I_exact_T_single_δ = @showprogress pmap(
    (K, T) -> F_I_exact_molecule(single_δ, 1, K, T),
    Ks_exact_mat,
    Ts_exact_mat,
)
writedlm("data/molecule/F_I_exact_T_single_delta.dat", res_F_I_exact_T_single_δ)

res_E_I_exact_T_single_δ = @showprogress pmap(
    (K, T) ->
        F_I_exact_molecule(single_δ, 1, K, T) -
        neg_TS_I_molecule(single_δ, 1, K, T),
    Ks_exact_mat,
    Ts_exact_mat,
)
writedlm("data/molecule/E_I_exact_T_single_delta.dat", res_E_I_exact_T_single_δ)

# QFT
res_F_I_QFT_T_single_M = @showprogress pmap(
    (K, T) -> F_I_molecule(0, single_M, K, T),
    Ks_QFT_mat,
    Ts_QFT_mat,
)
writedlm("data/molecule/F_I_QFT_T_single_M.dat", res_F_I_QFT_T_single_M)

res_F_I_QFT_T_single_δ = @showprogress pmap(
    (K, T) -> F_I_molecule(single_δ, 1, K, T),
    Ks_QFT_mat,
    Ts_QFT_mat,
)
writedlm("data/molecule/F_I_QFT_T_single_delta.dat", res_F_I_QFT_T_single_δ)

res_E_I_QFT_T_single_δ = @showprogress pmap(
    (K, T) -> E_I_molecule(single_δ, 1, K, T),
    Ks_QFT_mat,
    Ts_QFT_mat,
)
writedlm("data/molecule/E_I_QFT_T_single_delta.dat", res_E_I_QFT_T_single_δ)
