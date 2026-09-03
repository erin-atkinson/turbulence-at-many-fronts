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