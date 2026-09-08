# Just for this, background PV flux in gauge with total Jy=0
@inline function pv_flux_background_func(i, j, k, grid, sp, mean_fields, U)
    u = mean_fields.u
    v = mean_fields.v
    w = mean_fields.w
    b = mean_fields.b

    qxᶜᶠᶠ = ∂xᶜᶠᶠ(i, j, k, grid, q_func, sp, u, v, w, b)
    Uᶜᶠᶠ = ℑxyzᶜᶠᶠ(i, j, k, grid, U)

    return Uᶜᶠᶠ * qxᶜᶠᶠ
end

@doc raw"""
    PVFluxBackground(mean_fields, U, sp)
Return a kernel function operation that calculates the potential vorticity flux due to the background flow in the gauge with no along-front PV flux

```math
J_b = Uq - \int_{-\infty}^x \text{d}x\,U_x(x', z, t)q(x', z, t)  = \int_{-\infty}^x \text{d}x\,U(x', z, t)\frac{\partial q(x', z, t)}{\partial x}
```

See also [`PotentialVorticity`](@ref)
"""
function PVFluxBackground(mean_fields, U, sp)
    grid = mean_fields.u.grid
    loc = locationornothing((Face, Center, Face), mean_fields.u)
    return CumulativeIntegral(KernelFunctionOperation{loc...}(pv_flux_background_func, grid, sp, mean_fields, U); dims=1)
end

@inline function ref_to_zero(i, j, k, grid, field, loc)
    x, z = node(i, j, k, grid, loc...)
    
    return @inbounds field[i, j, k] - Oceananigans.Fields.interpolate((zero(x), z), field, loc, grid)
end

function RefToZero(field)
    grid = field.grid
    loc = location(field)
    return KernelFunctionOperation{loc...}(ref_to_zero, grid, field, (l() for l in loc))
end