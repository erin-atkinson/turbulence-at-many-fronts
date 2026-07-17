@inline function dsp_density_func(i, j, k, grid, clock, velocities, velocities_next, sp)
    t = clock.time

    u = velocities.u
    v = velocities.v
    
    u_next = velocities_next.u
    v_next = velocities_next.v

    αuu = ℑxᶜᵃᵃ(i, j, k, grid, αff_avg, (Face(), Center(), Center()), t, u, u_next, sp)
    αvv = ℑyᵃᶜᵃ(i, j, k, grid, αff_avg, (Center(), Face(), Center()), t, v, v_next, sp)
    
    return αuu - αvv
end

"""
    DSPDensity(u, v, u_next, v_next, U)
Deformation shear production density due to a background U
"""
function DSPDensity(clock, velocities, velocities_next, sp)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(dsp_density_func, grid, clock, velocities, velocities_next, sp)
end
