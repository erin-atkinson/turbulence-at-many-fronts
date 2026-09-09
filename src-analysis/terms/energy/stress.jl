@inline function stress_func(i, j, k, grid, clock, velocities, velocities_prev, sp)
    k = grid.Nz
    
    u = velocities.u
    v = velocities.v
    
    u_prev = velocities_prev.u
    v_prev = velocities_prev.v

    τu = ℑxᶜᵃᵃ(i, j, k, grid, f_avg_Gg, u, u_prev, u_flux_func, clock, sp)
    τv = ℑyᵃᶜᵃ(i, j, k, grid, f_avg_Gg, v, v_prev, v_flux_func, clock, sp)

    return -(τu + τv)
end

@doc raw"""
    STRESS(clock, velocities, velocities_prev, sp)
Return a kernel function operation that calculates work done by the velocity boundary conditions on the mean kinetic energy

The wind forcing is calcuated as
```math
\text{WIND} = \int \,\text{d}x \left [- \overline {\vec \tau}\cdot{\color{red} \overline{\vec u}(z=0)}\right]
```

This is a component of the mean kinetic energy equation:
```math
\frac{\text{d}}{\text{d}t}\text{MKE} = \text{DSP} + \text{WIND} + \text{BUOYANCY} + \text{SPONGE}_\text{MKE} - \text{LSP} - \text{VSP} + \text{STRAIN}_\text{MKE}
```

See also [`MKEDensity`](@ref), [`DSPDensity`](@ref), [`STRESS`](@ref), [`BUOYANCYDensity`](@ref), [`SPONGEMKEDensity`](@ref), [`LSPDensity`](@ref), [`VSPDensity`](@ref), [`STRAINMKEDensity`](@ref)
"""
function STRESS(clock, velocities, velocities_prev, sp)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Nothing), velocities.u)
    return KernelFunctionOperation{loc...}(stress_func, grid, clock, velocities, velocities_prev, sp)
end
