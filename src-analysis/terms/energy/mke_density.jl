@inline function mke_density_func(i, j, k, grid, velocities)
    u = velocities.u
    v = velocities.v
    w = velocities.w

    uu = ℑxᶜᵃᵃ(i, j, k, grid, fg, u, u)
    vv = ℑyᵃᶜᵃ(i, j, k, grid, fg, v, v)
    ww = ℑzᵃᵃᶜ(i, j, k, grid, fg, w, w)

    return (uu + vv + ww) / 2
end

"""
    MKEDensity(u, v, w)
Mean kinetic energy
"""
function MKEDensity(velocities)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(mke_density_func, grid, velocities)
end
