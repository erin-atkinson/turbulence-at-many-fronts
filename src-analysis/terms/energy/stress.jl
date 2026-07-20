@inline function τu_func(i, j, k, grid, loc, t, u, u_next, sp)
    x, y, z = node(i, j, k, grid, loc...)

    u_avg  = a_avg(i, j, k, grid, u, u_next)
    τx = -u_flux_func(x, t, sp) 

    return @inbounds τx * u_avg
end

@inline function τv_func(i, j, k, grid, loc, t, v, v_next, sp)
    x, y, z = node(i, j, k, grid, loc...)

    v_avg  = a_avg(i, j, k, grid, v, v_next)
    τy = -v_flux_func(x, t, sp) 

    return @inbounds τy * v_avg
end


@inline function stress_func(i, j, k, grid, clock, velocities, velocities_next, sp)
    t = clock.time
    k = grid.Nz
    
    u = velocities.u
    v = velocities.v
    
    u_next = velocities_next.u
    v_next = velocities_next.v

    τu = ℑxᶜᵃᵃ(i, j, k, grid, τu_func, (Face(), Center(), Center()), t, u, u_next, sp)
    τv = ℑyᵃᶜᵃ(i, j, k, grid, τv_func, (Center(), Face(), Center()), t, v, v_next, sp)

    return τu + τv
end

"""
    STRESS(clock, velocities, velocities_next, sp)
Wind stress forcing
"""
function STRESS(clock, velocities, velocities_next, sp)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Nothing), velocities.u)
    return KernelFunctionOperation{loc...}(stress_func, grid, clock, velocities, velocities_next, sp)
end