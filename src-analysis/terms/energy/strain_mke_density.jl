@inline function αff_avg(i, j, k, grid, loc, t, u, u_next, sp)
    x, y, z = node(i, j, k, grid, loc...)

    α = variable_strain_rate(t, sp) * strain_profile(x, sp)
    u_avg  = a_avg(i, j, k, grid, u, u_next)

    return @inbounds -α * u[i, j, k] * u_avg
end


@inline function strain_mke_density_func(i, j, k, grid, clock, velocities, velocities_next, sp)
    t = clock.time

    u = velocities.u
    v = velocities.v
    w = velocities.w
    
    u_next = velocities_next.u
    v_next = velocities_next.v
    w_next = velocities_next.w

    αuu = ℑxᶜᵃᵃ(i, j, k, grid, αff_avg, (Face(), Center(), Center()), t, u, u_next, sp)
    αvv = ℑyᵃᶜᵃ(i, j, k, grid, αff_avg, (Center(), Face(), Center()), t, v, v_next, sp)
    αww = ℑzᵃᵃᶜ(i, j, k, grid, αff_avg, (Center(), Center(), Face()), t, w, w_next, sp)
    
    return (αuu + αvv + αww) / 2
end

"""
    STRAINMKEDensity(velocities, velocities_next, turbulent_fluxes)
Input coarse-grained velocities and fluxes for lateral shear production density
"""
function STRAINMKEDensity(clock, velocities, velocities_next, sp)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(strain_mke_density_func, grid, clock, velocities, velocities_next, sp)
end
