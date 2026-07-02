# Coarse-grained kernel sizes
coarse_σx = sp.Lh / 75
coarse_σz = 4.0

const weno = WENO(; order=9, weight_computation=Oceananigans.Utils.BackendOptimizedDivision)
const centered = Centered(; order=2)
