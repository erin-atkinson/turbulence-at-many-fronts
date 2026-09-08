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
Return a kernel function operation that calculates the work done by the strain flow on the resolved flow

The deformation shear production is calcuated as
```math
\text{DSP} = \int \,\text{d}x\text{d}z \left [ - \frac{\partial U}{\partial x}(\overline{u}{\color{red} \overline{u}} - \overline{v}{\color{red} \overline{v}})\right]
```

This is a component of the mean kinetic energy equation:
```math
\frac{\text{d}}{\text{d}t}\text{MKE} = \text{DSP} + \text{WIND} + \text{BUOYANCY} + \text{SPONGE}_\text{MKE} - \text{LSP} - \text{VSP} + \text{STRAIN}_\text{MKE}
```

See also [`MKEDensity`](@ref), [`STRESS`](@ref), [`BUOYANCYDensity`](@ref), [`SPONGEMKEDensity`](@ref), [`LSPDensity`](@ref), [`VSPDensity`](@ref), [`STRAINMKEDensity`](@ref)
"""
function DSPDensity(clock, velocities, velocities_prev, sp)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(dsp_density_func, grid, clock, velocities, velocities_prev, sp)
end
