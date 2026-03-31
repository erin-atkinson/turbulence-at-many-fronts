@inline function KE_func(i, j, k, grid, clock, u, v, w)
    
    U = ℑxᶜᵃᵃ(i, j, k, grid, u)
    V = ℑyᵃᶜᵃ(i, j, k, grid, v)
    W = ℑzᵃᵃᶜ(i, j, k, grid, w)

    KE = (U^2 + V^2 + W^2) / 2
    return KE
end

"""
    KE(u, v, w)
Kinetic energy
"""
function KE(u, v, w)
    grid = u.grid
    loc = locationornothing((Center, Center, Center), u)
    return KernelFunctionOperation{loc...}(KE_func, grid, u, v, w)
end
