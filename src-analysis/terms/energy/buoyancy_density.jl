@inline function buoyancy_density_func(i, j, k, grid, w, w_prev, b)
    return @inbounds ℑzᵃᵃᶜ(i, j, k, grid, a_avg, w, w_prev) * b[i, j, k]
end

@doc raw"""
    BUOYANCYDensity(w, w_prev, b)
Production of mean kinetic energy from potential energy of the mean state
```math
\text{BUOYANCY} = \int \text{d}V \overline{w}\overline{b}
```
"""
function BUOYANCYDensity(w, w_prev, b)
    grid = w.grid
    loc = locationornothing((Center, Center, Center), w)
    return KernelFunctionOperation{loc...}(buoyancy_density_func, grid, w, w_prev, b)
end
