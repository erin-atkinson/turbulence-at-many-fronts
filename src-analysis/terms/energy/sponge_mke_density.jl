@inline function σfg(i, j, k, grid, loc, u, u_next, sp)
    x, y, z = node(i, j, k, grid, loc...)

    u_avg  = a_avg(i, j, k, grid, u, u_next)

    σ = sponge_layer_func(z, sp)

    return @inbounds -σ * u[i, j, k] * u_avg
end


@inline function sponge_mke_density_func(i, j, k, grid, velocities, velocities_next, sp)
    u = velocities.u
    v = velocities.v
    w = velocities.w
    
    u_next = velocities_next.u
    v_next = velocities_next.v
    w_next = velocities_next.w

    σuu = ℑxᶜᵃᵃ(i, j, k, grid, σfg, (Face(), Center(), Center()), u, u_next, sp)
    σvv = ℑyᵃᶜᵃ(i, j, k, grid, σfg, (Center(), Face(), Center()), v, v_next, sp)
    σww = ℑzᵃᵃᶜ(i, j, k, grid, σfg, (Center(), Center(), Face()), w, w_next, sp)

    return σuu + σvv + σww
end

"""
    SPONGEMKEDensity(velocities, velocities_next, sp)
Sponge layer forcing
"""
function SPONGEMKEDensity(velocities, velocities_next, sp)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(sponge_mke_density_func, grid, velocities, velocities_next, sp)
end