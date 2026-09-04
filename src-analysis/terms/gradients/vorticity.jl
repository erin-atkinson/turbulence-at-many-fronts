@inline function vorticity_x_func(i, j, k, grid, v, w)
    return ∂yᶜᶠᶠ(i, j, k, grid, w) - ∂zᶜᶠᶠ(i, j, k, grid, v)
end

@inline function vorticity_y_func(i, j, k, grid, u, w)
    return ∂zᶠᶜᶠ(i, j, k, grid, u) - ∂xᶠᶜᶠ(i, j, k, grid, w)
end

@inline function vorticity_z_func(i, j, k, grid, u, v)
    return ∂xᶠᶠᶜ(i, j, k, grid, v) - ∂yᶠᶠᶜ(i, j, k, grid, u)
end

@doc raw"""
    VorticityX(v, w)
Return a kernel function operation that calculates the vorticity in the x direction.

See also [`VorticityY`](@ref), [`VorticityZ`](@ref), [`PotentialVorticity`](@ref)

```math
\omega_x = \frac{\partial w}{\partial y} - \frac{\partial v}{\partial z}
```
"""
function VorticityX(v, w)
    grid = v.grid
    loc = locationornothing((Center, Face, Face), v)
    return KernelFunctionOperation{loc...}(vorticity_x_func, grid, v, w)
end


@doc raw"""
    VorticityY(v, w)
Return a kernel function operation that calculates the vorticity in the y direction.

See also [`VorticityX`](@ref), [`VorticityZ`](@ref), [`PotentialVorticity`](@ref)

```math
\omega_y = \frac{\partial u}{\partial z} - \frac{\partial w}{\partial x}
```
"""
function VorticityY(u, w)
    grid = u.grid
    loc = locationornothing((Face, Center, Face), u)
    return KernelFunctionOperation{loc...}(vorticity_y_func, grid, u, w)
end


@doc raw"""
    VorticityZ(v, w)
Return a kernel function operation that calculates the vorticity in the z direction.

See also [`VorticityX`](@ref), [`VorticityY`](@ref), [`PotentialVorticity`](@ref)

```math
\omega_z = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}
```
"""
function VorticityZ(u, v)
    grid = v.grid
    loc = locationornothing((Face, Face, Center), v)
    return KernelFunctionOperation{loc...}(vorticity_z_func, grid, u, v)
end

VorticityX(_, v, w) = VorticityX(v, w)
VorticityY(u, _, w) = VorticityY(u, w)
VorticityZ(u, v, _) = VorticityZ(u, v)
