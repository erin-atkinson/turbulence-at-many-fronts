@inline function bflux_density_func(i, j, k, grid, wb)
    return ℑzᵃᵃᶜ(i, j, k, grid, wb)
end

@doc raw"""
    BFLUXDensity(wb)
Return a kernel function operation that calculates the production of turbulent kinetic energy from mean potential energy.

The turbulent buoyancy flux is calculated as
```math
\text{BFLUX} = \int \text{d}V w'b'
```

This is a component of the mean potential energy equation:
```math
\frac{\text{d}}{\text{d}t}\text{MPE} = -\text{BUOYANCY} + \text{SPONGE}_\text{MPE} - \text{BFLUX} + \text{MIXED} + \text{COOLING} + \text{STRAIN}_\text{MPE}
```

See also [`MPEDensity`](@ref), [`MLD`](@ref), [`BUOYANCYDensity`](@ref), [`SPONGEMPEDensity`](@ref), [`MIXEDDensity`](@ref), [`COOLING`](@ref), [`STRAINMPEDensity`](@ref)
"""
function BFLUXDensity(wb)
    grid = wb.grid
    loc = locationornothing((Center, Center, Center), wb)
    return KernelFunctionOperation{loc...}(bflux_density_func, grid, wb)
end
