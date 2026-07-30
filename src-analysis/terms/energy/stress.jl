@inline function τu_func(i, j, k, grid, loc, t, u, u_prev, sp)
    x, y, z = node(i, j, k, grid, loc...)

    u_avg  = a_avg(i, j, k, grid, u, u_prev)
    τx = -u_flux_func(x, t, sp) 

    return @inbounds τx * u_avg
end

@inline function τv_func(i, j, k, grid, loc, t, v, v_prev, sp)
    x, y, z = node(i, j, k, grid, loc...)

    v_avg  = a_avg(i, j, k, grid, v, v_prev)
    τy = -v_flux_func(x, t, sp) 

    return @inbounds τy * v_avg
end


@inline function stress_func(i, j, k, grid, clock, velocities, velocities_prev, sp)
    t = clock.time
    k = grid.Nz
    
    u = velocities.u
    v = velocities.v
    
    u_prev = velocities_prev.u
    v_prev = velocities_prev.v

    τu = ℑxᶜᵃᵃ(i, j, k, grid, τu_func, (Face(), Center(), Center()), t, u, u_prev, sp)
    τv = ℑyᵃᶜᵃ(i, j, k, grid, τv_func, (Center(), Face(), Center()), t, v, v_prev, sp)

    return τu + τv
end

"""
    STRESS(clock, velocities, velocities_prev, sp)
Wind stress forcing
"""
function STRESS(clock, velocities, velocities_prev, sp)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Nothing), velocities.u)
    return KernelFunctionOperation{loc...}(stress_func, grid, clock, velocities, velocities_prev, sp)
end