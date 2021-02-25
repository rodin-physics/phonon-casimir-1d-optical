include("../src/molecule.jl")

## Parameters
# Masses are given in units of m, the atomic mass
# Potentials are given in units of k
# Temperatres are given in units of Ω = √(k / m)

nPts_exact = 2000;  # Number of points used in Exact Diagonalization plots
nPts_QFT = 20;      # Number of points used in QFT plots

M_min = 1 / 20;     # Minimum impurity mass
M_max = 10000;      # Maximum impurity mass

δ_min = 0;          # Minimum external potential
δ_max = 100;        # Maximum external potential

single_M = 10;      # Impurity mass values used in showing T dependence
single_δ = 1;       # External potential value used in showing T dependence

T_min = 0;          # Minimum temperature
T_max = 1;          # Maximum temperature

Ks = [0.001, 0.01, 0.1, 1, 10]; # Confining potential values

## Computed parameters

# Exact diagonalization parameters
Ms_exact = (range(M_min, M_max, length = nPts_exact))
# NOTE: log_Ms is a log-spaced array of Ms, not ln(Ms)
log_Ms_exact = exp.(range(log.(M_min), log.(M_max), length = nPts_exact))
δs_exact = range(δ_min, δ_max, length = nPts_exact)
Ts_exact = range(T_min, T_max, length = nPts_exact)

# QFT parameters
Ms_QFT = (range(M_min, M_max, length = nPts_QFT))
# NOTE: log_Ms is a log-spaced array of Ms, not ln(Ms)
log_Ms_QFT = exp.(range(log.(M_min), log.(M_max), length = nPts_QFT))
δs_QFT = range(δ_min, δ_max, length = nPts_QFT)
Ts_QFT = range(T_min, T_max, length = nPts_QFT)

# Matrices used in pmap
Ms_exact_mat = repeat(Ms_exact, 1, length(Ks))
log_Ms_exact_mat = repeat(log_Ms_exact, 1, length(Ks))
δs_exact_mat = repeat(δs_exact, 1, length(Ks))
Ts_exact_mat = repeat(Ts_exact, 1, length(Ks))
Ks_exact_mat = repeat(Ks', nPts_exact, 1)

Ms_QFT_mat = repeat(Ms_QFT, 1, length(Ks))
log_Ms_QFT_mat = repeat(log_Ms_QFT, 1, length(Ks))
δs_QFT_mat = repeat(δs_QFT, 1, length(Ks))
Ts_QFT_mat = repeat(Ts_QFT, 1, length(Ks))
Ks_QFT_mat = repeat(Ks', nPts_QFT, 1)
