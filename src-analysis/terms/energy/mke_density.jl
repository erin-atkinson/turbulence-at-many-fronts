@inline function mke_density_func(i, j, k, grid, u, v, w)
    
    uu = ℑxᶜᵃᵃ(i, j, k, grid, fg, u, u)
    vv = ℑyᵃᶜᵃ(i, j, k, grid, fg, v, v)
    ww = ℑzᵃᵃᶜ(i, j, k, grid, fg, w, w)

    return (uu + vv + ww) / 2
end

"""
    MKE(u, v, w)
Mean kinetic energy
"""
function MKE(u, v, w)
    grid = u.grid
    loc = locationornothing((Center, Center, Center), u)
    return KernelFunctionOperation{loc...}(mke_density_func, grid, u, v, w)
end
