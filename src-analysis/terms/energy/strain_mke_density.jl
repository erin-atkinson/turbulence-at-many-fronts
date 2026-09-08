@inline function strain_mke_density_func(i, j, k, grid, clock, velocities, velocities_prev, sp)
    t = clock.time

    u = velocities.u
    v = velocities.v
    w = velocities.w
    
    u_prev = velocities_prev.u
    v_prev = velocities_prev.v
    w_prev = velocities_prev.w

    αuu = ℑxᶜᵃᵃ(i, j, k, grid, αff_avg, (Face(), Center(), Center()), t, u, u_prev, sp)
    αvv = ℑyᵃᶜᵃ(i, j, k, grid, αff_avg, (Center(), Face(), Center()), t, v, v_prev, sp)
    αww = ℑzᵃᵃᶜ(i, j, k, grid, αff_avg, (Center(), Center(), Face()), t, w, w_prev, sp)
    
    return -(αuu + αvv + αww) / 2
end

@doc raw"""
    STRAINMKEDensity(clock, velocities, velocities_prev, sp)
Return a kernel function operation that calculates the change in mean kinetic energy due to the strain flow.

Note that this is not the deformation shear production. This term represents the net flux of kinetic energy into the domain by the strain flow. This is calcuated as
```math
\text{STRAIN}_\text{MKE} = \int \,\text{d}x\text{d}z \left [ \frac{\partial U}{\partial x}\frac{\overline{\vec{u}}\cdot {\color{red} \overline{\vec{u}}}}{2}\right]
```

This is a component of the mean kinetic energy equation:
```math
\frac{\text{d}}{\text{d}t}\text{MKE} = \text{DSP} + \text{WIND} + \text{BUOYANCY} + \text{SPONGE}_\text{MKE} - \text{LSP} - \text{VSP} + \text{STRAIN}_\text{MKE}
```

See also [`MKEDensity`](@ref), [`DSPDensity`](@ref), [`WINDDensity`](@ref), [`BUOYANCYDensity`](@ref), [`SPONGEMKEDensity`](@ref), [`LSPDensity`](@ref), [`VSPDensity`](@ref)
"""
function STRAINMKEDensity(clock, velocities, velocities_prev, sp)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(strain_mke_density_func, grid, clock, velocities, velocities_prev, sp)
end
