# Coarse-grained kernel sizes
coarse_σx = sp.Lh / 75
coarse_σz = sp.H / 25

const weno = WENO(; order=9)
const centered = Centered(; order=2)
