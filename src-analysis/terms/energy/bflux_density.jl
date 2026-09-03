@inline function bflux_density_func(i, j, k, grid, wb)
    return ℑzᵃᵃᶜ(i, j, k, grid, wb)
end

@doc raw"""
    BFLUXDensity(wb)
Production of turbulent kinetic energy from potential energy of the mean state
```math
\text{BFLUX} = \int \text{d}V w'b'
```
"""
function BFLUXDensity(wb)
    grid = wb.grid
    loc = locationornothing((Center, Center, Center), wb)
    return KernelFunctionOperation{loc...}(bflux_density_func, grid, wb)
end
