@inline function TKE_func(i, j, k, grid, clock, u, v, w, U, V, W)

    u = ℑxᶜᵃᵃ(i, j, k, grid, u)
    v = ℑyᵃᶜᵃ(i, j, k, grid, v)
    w = ℑzᵃᵃᶜ(i, j, k, grid, w)
    
    U = ℑxᶜᵃᵃ(i, j, k, grid, U)
    V = ℑyᵃᶜᵃ(i, j, k, grid, V)
    W = ℑzᵃᵃᶜ(i, j, k, grid, V)

    KE = (u^2 + v^2 + w^2) / 2
    MKE = (U^2 + V^2 + W^2) / 2
    return KE - MKE
end

"""
    TKE(u, v, w)
Turbulent kinetic energy
"""
function KE(u, v, w, U, V, W)
    grid = u.grid
    loc = locationornothing((Center, Center, Center), u)
    return KernelFunctionOperation{loc...}(KE_func, grid, u, v, w)
end
