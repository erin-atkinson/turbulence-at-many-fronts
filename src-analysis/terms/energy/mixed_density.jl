@inline function mixed_density_func(i, j, k, grid, clock, b, b_prev, h_ml, h_ml_prev)
    Δt = clock.last_Δt

    b_avg = a_avg(i, j, k, grid, b, b_prev)
    ∂h_ml∂t = @inbounds (h_ml[i, j, k] - h_ml_prev[i, j, k]) / Δt

    return -b_avg * ∂h_ml∂t
end

@doc raw"""
    MIXEDDensity(clock, b, b_prev, h_ml, h_ml_prev)
Return a kernel function operation that calculates the change of potential energy of the mixed layer due to a change in mixed layer depth

This is a component of the mean potential energy equation:
```math
\frac{\text{d}}{\text{d}t}\text{MPE} = -\text{BUOYANCY} + \text{SPONGE}_\text{MPE} - \text{BFLUX} + \text{MIXED} + \text{COOLING} + \text{STRAIN}_\text{MPE}
```

See also [`MPEDensity`](@ref), [`MLD`](@ref), [`BUOYANCYDensity`](@ref), [`SPONGEMPEDensity`](@ref), [`BFLUXDensity`](@ref), [`COOLING`](@ref), [`STRAINMPEDensity`](@ref)
"""
function MIXEDDensity(clock, b, b_prev, h_ml, h_ml_prev)
    grid = b.grid
    loc = locationornothing((Center, Center, Center), b)
    return KernelFunctionOperation{loc...}(mixed_density_func, grid, clock, b, b_prev, h_ml, h_ml_prev)
end
