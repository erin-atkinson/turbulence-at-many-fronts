@inline function buoyancy_density_func(i, j, k, grid, w, w_prev, b)
    return @inbounds ℑzᵃᵃᶜ(i, j, k, grid, a_avg, w, w_prev) * b[i, j, k]
end

@doc raw"""
    BUOYANCYDensity(w, w_prev, b)
Return a kernel function operation that calculates the production of mean kinetic energy from mean potential energy.

The buoyancy flux is calculated as
```math
\text{BFLUX} = \int \text{d}V \overline{w}\overline{b}
```

This is a component of the mean kinetic energy equation:
```math
\frac{\text{d}}{\text{d}t}\text{MKE} = \text{DSP} + \text{WIND} + \text{BUOYANCY} + \text{SPONGE}_\text{MKE} - \text{LSP} - \text{VSP} + \text{STRAIN}_\text{MKE}
```

as well as the mean potential energy equation:
```math
\frac{\text{d}}{\text{d}t}\text{MPE} = -\text{BUOYANCY} + \text{SPONGE}_\text{MPE} - \text{BFLUX} + \text{MIXED} + \text{COOLING} + \text{STRAIN}_\text{MPE}
```

See also [`MPEDensity`](@ref), [`MKEDensity`](@ref)
"""
function BUOYANCYDensity(w, w_prev, b)
    grid = w.grid
    loc = locationornothing((Center, Center, Center), w)
    return KernelFunctionOperation{loc...}(buoyancy_density_func, grid, w, w_prev, b)
end
