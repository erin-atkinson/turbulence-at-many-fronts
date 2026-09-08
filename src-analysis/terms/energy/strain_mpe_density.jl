@inline function strain_mpe_density_func(i, j, k, grid, clock, b, h_ml, h_ml_prev, sp)
    t = clock.time
    x, y, z = node(i, j, k, grid, Center(), Center(), Center())
    U = variable_strain_rate(t, sp) * velocity_profile(x, sp)

    h_ml_avg = a_avg(i, j, k, grid, h_ml, h_ml_prev)
    
    return U * ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b) * (z + h_ml_avg)
end

@doc raw"""
    STRAINMPEDensity(clock, b, h_ml, h_ml_prev, sp)
Return a kernel function operation that calculates the change in mean potential energy due to the strain flow

The result is calculated as follows:
```math
\text{STRAIN}_\text{MPE} = \int \,\text{d}x\text{d}z \left [ U\frac{\partial \overline{b}}{\partial x} (z +{\color{red} h_{ml}})\right]
```

This is a component of the mean potential energy equation:
```math
\frac{\text{d}}{\text{d}t}\text{MPE} = -\text{BUOYANCY} + \text{SPONGE}_\text{MPE} - \text{BFLUX} + \text{MIXED} + \text{COOLING} + \text{STRAIN}_\text{MPE}
```

See also [`MPEDensity`](@ref), [`MLD`](@ref), [`BUOYANCYDensity`](@ref), [`SPONGEMPEDensity`](@ref), [`BFLUXDensity`](@ref), [`MIXEDDensity`](@ref), [`COOLING`](@ref)
"""
function STRAINMPEDensity(clock, b, h_ml, h_ml_prev, sp)
    grid = b.grid
    loc = locationornothing((Center, Center, Center), b)
    return KernelFunctionOperation{loc...}(strain_mpe_density_func, grid, clock, b, h_ml, h_ml_prev, sp)
end
