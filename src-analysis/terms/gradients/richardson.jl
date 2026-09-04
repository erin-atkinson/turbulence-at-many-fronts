@inline function balanced_richardson_func(i, j, k, grid, sp, b)
    bx = ℑxzᶜᵃᶠ(i, j, k, grid, ∂xᶠᶜᶜ, b)
    by = ℑyzᵃᶜᶠ(i, j, k, grid, ∂yᶜᶠᶜ, b)
    bz = ∂zᶜᶜᶠ(i, j, k, grid, b)

    S² = (bx^2 + by^2) / sp.f^2

    return bz / S²
end

@doc raw"""
    BalancedRichardson(b, sp)
Return a kernel function operation that calculates the Richardson number of the balanced flow due to `b` only.

See also 
[`Richardson`](@ref)

```math
\text{Ri}_b = \frac{N^2}{S^2}
```
where 
```math
N^2 = \frac{\partial b}{\partial z} \quad \text{and} \quad S^2 = \frac{1}{f^2}\left (\frac{\partial b}{\partial x}^2 + \frac{\partial b}{\partial x}^2 \right)
```
"""
function BalancedRichardson(b, sp)
    grid = b.grid
    loc = locationornothing((Center, Center, Face), b)
    return KernelFunctionOperation{loc...}(balanced_richardson_func, grid, sp, b)
end

@inline function richardson_func(i, j, k, grid, u, v, b)
    Sx = ℑxᶜᵃᵃ(i, j, k, grid, ∂zᶠᶜᶠ, u)
    Sy = ℑyᵃᶜᵃ(i, j, k, grid, ∂zᶜᶠᶠ, v)
    bz = ∂zᶜᶜᶠ(i, j, k, grid, b)
    
    return bz / (Sx^2 + Sy^2)
end

@doc raw"""
    Richardson(u, v, b)
Return a kernel function operation that calculates the Richardson number of the flow.

See also 
[`BalancedRichardson`](@ref)

```math
\text{Ri} = \frac{N^2}{S^2}
```
where 
```math
N^2 = \frac{\partial b}{\partial z} \quad \text{and} \quad S^2 = \frac{\partial u}{\partial z}^2 + \frac{\partial v}{\partial z}^2
```
```
"""
function Richardson(u, v, b)
    grid = u.grid
    loc = locationornothing((Center, Center, Face), u)
    return KernelFunctionOperation{loc...}(richardson_func, grid, u, v, b)
end
