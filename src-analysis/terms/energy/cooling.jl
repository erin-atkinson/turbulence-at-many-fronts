function cooling_func(i, j, k, grid, clock, h_ml, h_ml_prev, sp)
    t = clock.time
    B = b_flux_func(t, sp)
    return @inbounds a_avg(i, j, k, grid, h_ml, h_ml_prev) * B * sp.Lx 
end

@doc raw"""
    Cooling(clock, h_ml, h_ml_prev, sp)
Return a kernel function operation that calculates potential energy input by the surface cooling, referenced to the base of the mixed layer

This is a component of the mean potential energy equation:
```math
\frac{\text{d}}{\text{d}t}\text{MPE} = -\text{BUOYANCY} + \text{SPONGE}_\text{MPE} - \text{BFLUX} + \text{MIXED} + \text{COOLING} + \text{STRAIN}_\text{MPE}
```

See also [`MPEDensity`](@ref), [`MLD`](@ref), [`BUOYANCYDensity`](@ref), [`SPONGEMPEDensity`](@ref), [`BFLUXDensity`](@ref), [`MIXEDDensity`](@ref), [`STRAINMPEDensity`](@ref)
"""
function COOLING(clock, h_ml, h_ml_prev, sp)
    grid = h_ml.grid
    return KernelFunctionOperation{Nothing, Nothing, Nothing}(cooling_func, grid, clock, h_ml, h_ml_prev, sp)
end
