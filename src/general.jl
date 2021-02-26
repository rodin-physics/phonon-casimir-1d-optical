using DelimitedFiles
using LaTeXStrings
using LinearAlgebra
using Plots
using ProgressMeter
using QuadGK

## Parameters
const ν = 1e-5;             # Relative tolerance for integration
const η = 1e-5;             # Small number used for i0
const α = 1e-9;             # Small number for absolute tolerance
const max_omega = 10000;    # Large number used to truncate the finite
# temperature Matsubara sum
const NumEvals = 1e7;       # Maximum number of evaluations in quadgk

# Colors for plotting
my_red = RGB(215 / 255, 67 / 255, 84 / 255)
my_green = RGB(106 / 255, 178 / 255, 71 / 255)
my_blue = RGB(100 / 255, 101 / 255, 218 / 255)
my_violet = RGB(169 / 255, 89 / 255, 201 / 255)
my_orange = RGB(209 / 255, 135 / 255, 46 / 255)

colors = [my_red, my_green, my_blue, my_violet, my_orange]
