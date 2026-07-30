@inline function σfg(i, j, k, grid, loc, u, u_prev, sp)
    x, y, z = node(i, j, k, grid, loc...)

    u_avg  = a_avg(i, j, k, grid, u, u_prev)

    σ = sponge_layer_func(z, sp)

    return @inbounds -σ * u[i, j, k] * u_avg
end


@inline function sponge_mke_density_func(i, j, k, grid, velocities, velocities_prev, sp)
    u = velocities.u
    v = velocities.v
    w = velocities.w
    
    u_prev = velocities_prev.u
    v_prev = velocities_prev.v
    w_prev = velocities_prev.w

    σuu = ℑxᶜᵃᵃ(i, j, k, grid, σfg, (Face(), Center(), Center()), u, u_prev, sp)
    σvv = ℑyᵃᶜᵃ(i, j, k, grid, σfg, (Center(), Face(), Center()), v, v_prev, sp)
    σww = ℑzᵃᵃᶜ(i, j, k, grid, σfg, (Center(), Center(), Face()), w, w_prev, sp)

    return σuu + σvv + σww
end

"""
    SPONGEMKEDensity(velocities, velocities_prev, sp)
Sponge layer forcing
"""
function SPONGEMKEDensity(velocities, velocities_prev, sp)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(sponge_mke_density_func, grid, velocities, velocities_prev, sp)
end