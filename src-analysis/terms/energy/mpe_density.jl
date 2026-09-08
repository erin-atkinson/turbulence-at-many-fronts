@inline function mpe_density_func(i, j, k, grid, b, h_ml)
    x, y, z = node(i, j, k, grid, Center(), Center(), Center())
    
    return @inbounds -b[i, j, k] * (z + h_ml[i, j, k])
end

@doc raw"""
    MPEDensity(b, h_ml)
Return a kernel function operation that calculates the mean potential energy referenced to the base of the mixed layer

The mean potential energy is calculated as follows:
```math
\text{MPE} = -\int\text{d}x\text{d}z\, \overline{b}(z + h)
```

This is a component of the mean potential energy equation:
```math
\frac{\text{d}}{\text{d}t}\text{MPE} = -\text{BUOYANCY} + \text{SPONGE}_\text{MPE} - \text{BFLUX} + \text{MIXED} + \text{COOLING} + \text{STRAIN}_\text{MPE}
```

See also [`MLD`](@ref), [`BUOYANCYDensity`](@ref), [`SPONGEMPEDensity`](@ref), [`BFLUXDensity`](@ref), [`MIXEDDensity`](@ref), [`COOLING`](@ref), [`STRAINMPEDensity`](@ref)
"""
function MPEDensity(b, h_ml)
    grid = b.grid
    loc = locationornothing((Center, Center, Center), b)
    return KernelFunctionOperation{loc...}(mpe_density_func, grid, b, h_ml)
end
