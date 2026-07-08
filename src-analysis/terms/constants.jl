# Coarse-grained kernel sizes
coarse_σx = 10.0
coarse_σz = 1.0

const weno = WENO(; order=9, weight_computation=Oceananigans.Utils.BackendOptimizedDivision)
const centered = Centered(; order=2)
