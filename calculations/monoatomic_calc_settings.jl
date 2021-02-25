include("../src/chain.jl")
## Parameters
ds = 1:1:20         # Separations between defects
N = 1000;           # Number of unit cells in the chain (each contains 2 atoms)
T = 0.02;           # Temperature for finite-T plots
K = 1e-6;           # Confining potential
m = 1;              # Mass in units of μ
Ms = [m, m]         # Array of masses for the diatomic chain

# Defect parameters
Defect_Ms = [0.1, 0.5, 5, 50];  # Impurity masses used in mass-mass interaction
Defect_Δs = [0.1, 0.2, 0.5, 5]; # Δs used in well-well interactions
Imp_M = 10;     # Impurity mass used in mixed-defect calculations
Imp_Δ = 10      # Δ used in mixed-defect calculations

# Cluster parameters
left_imp = Impurity(1, 1, m, 5)
right_imp = Impurity(10, 1, m, 5)
Ms_cluster = [0.5, 0.8, 1.5, 2.0]
# Distances of the middle impurity from the left one
Ds = 9:14

# Heatmap parameters
nPts_heatmap = 100; # Number of points along each direction in the heatmap
MΔ_max = 5;         # Maximum value of M and Δ for the heatmap


## Computed Parameters

MΔ = range(1e-5, MΔ_max, length = nPts_heatmap)
MΔ_mat = repeat(MΔ, 1, nPts_heatmap)

# Matrices for pmap
Ds_cluster_mat = repeat(Ds, 1, length(Ms_cluster))
Ms_cluster_mat = repeat(Ms_cluster', length(Ds), 1)

Defect_Ms_mat = repeat(Defect_Ms, 1, length(ds))
Defect_Δs_mat = repeat(Defect_Δs, 1, length(ds))

ds_M_mat = repeat(ds',length(Defect_Ms), 1)
ds_Δ_mat = repeat(ds',length(Defect_Δs), 1)
