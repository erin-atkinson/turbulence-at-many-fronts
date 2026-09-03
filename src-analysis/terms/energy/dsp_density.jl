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

@doc raw"""
    DSPDensity(clock, velocities, velocities_prev, sp)
Deformation shear production density due to the strain flow
```math
\text{DSP} = \int \text{d}V \alpha (\overline{u}^2  - \overline{v}^2)
```
"""
function DSPDensity(clock, velocities, velocities_prev, sp)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(dsp_density_func, grid, clock, velocities, velocities_prev, sp)
end
