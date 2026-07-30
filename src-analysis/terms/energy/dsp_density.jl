@inline function dsp_density_func(i, j, k, grid, clock, velocities, velocities_prev, sp)
    t = clock.time

    u = velocities.u
    v = velocities.v
    
    u_prev = velocities_prev.u
    v_prev = velocities_prev.v

    αuu = ℑxᶜᵃᵃ(i, j, k, grid, αff_avg, (Face(), Center(), Center()), t, u, u_prev, sp)
    αvv = ℑyᵃᶜᵃ(i, j, k, grid, αff_avg, (Center(), Face(), Center()), t, v, v_prev, sp)
    
    return αuu - αvv
end

"""
    DSPDensity(u, v, u_prev, v_prev, U)
Deformation shear production density due to a background U
"""
function DSPDensity(clock, velocities, velocities_prev, sp)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(dsp_density_func, grid, clock, velocities, velocities_prev, sp)
end
