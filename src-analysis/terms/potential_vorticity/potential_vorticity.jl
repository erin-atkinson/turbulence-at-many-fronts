@inline function potential_vorticity_func(i, j, k, grid, sp, u, v, w, b)
    ηx = ℑxᶠᵃᵃ(i, j, k, grid, vorticity_x_func, v, w)
    ηy = ℑyᵃᶠᵃ(i, j, k, grid, vorticity_y_func, u, w)
    ηz = ℑzᵃᵃᶠ(i, j, k, grid, vorticity_z_func, u, v) + sp.f

    bx = ∂xᶠᶜᶜ(i, j, k, grid, b)
    by = ∂yᶜᶠᶜ(i, j, k, grid, b)
    bz = ∂zᶜᶜᶠ(i, j, k, grid, b)
    
    return (bx * ηx
          + by * ηy
          + bz * ηz)
end

@doc raw"""
    PotentialVorticity(u, v, w, b, sp)
Return a kernel function operation that calculates the potential vorticity.
```math
q = (\nabla \times \overline{\vec{u}} + f\hat z) \cdot \nabla \overline{b}
```
"""
function PotentialVorticity(u, v, w, b, sp)
    grid = u.grid
    loc = locationornothing((Face, Face, Face), u)
    return KernelFunctionOperation{loc...}(q_func, grid, sp, u, v, w, b)
end
