include("../src/chain.jl")
## Parameters
ds = 1:1:20         # Separations between defects
N = 1000;           # Number of unit cells in the chain (each contains 2 atoms)
T = 0.02;           # Temperature for finite-T plots
K = 1e-6;           # Confining potential
m1 = 1;             # Mass in units of μ
m2 = 3;             # Mass in units of μ
Ms = [m1, m2]       # Array of masses for the diatomic chain

Imp_M = 1.8;        # Impurity mass

# Cluster parameters
left_imp = Impurity(1, 1, Ms[1], 5)
right_imp = Impurity(10, 1, Ms[1], 5)
Ms_cluster = [0.5, 1.8, 4.2]
# Distances of the middle impurity from the left one
Ds = 9:14
# Matrices for pmap
Ds_cluster_mat = repeat(Ds, 1, length(Ms_cluster))
Ms_cluster_mat = repeat(Ms_cluster', length(Ds), 1)
