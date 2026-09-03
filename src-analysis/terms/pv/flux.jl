@inline function pv_flux_x_func(i, j, k, grid, sp, mean_fields)
    # x-flux of PV: cff
    
    u = mean_fields.u
    v = mean_fields.v
    w = mean_fields.w
    b = mean_fields.b

    qᶜᶠᶠ = ℑxᶜᵃᵃ(i, j, k, grid, q_func, sp, u, v, w, b)
    uᶜᶠᶠ = ℑxyzᶜᶠᶠ(i, j, k, grid, u)

    return qᶜᶠᶠ * uᶜᶠᶠ
end

@inline function pv_flux_y_func(i, j, k, grid, sp, mean_fields)
    # y-flux of PV: fcf
    
    u = mean_fields.u
    v = mean_fields.v
    w = mean_fields.w
    b = mean_fields.b

    qᶠᶜᶠ = ℑyᵃᶜᵃ(i, j, k, grid, q_func, sp, u, v, w, b)
    vᶠᶜᶠ = ℑxyzᶠᶜᶠ(i, j, k, grid, v)

    return qᶠᶜᶠ * vᶠᶜᶠ
end

@inline function pv_flux_z_func(i, j, k, grid, sp, mean_fields)
    # z-flux of PV: ffc
    
    u = mean_fields.u
    v = mean_fields.v
    w = mean_fields.w
    b = mean_fields.b

    qᶠᶠᶜ = ℑzᵃᵃᶜ(i, j, k, grid, q_func, sp, u, v, w, b)
    wᶠᶠᶜ = ℑxyzᶠᶠᶜ(i, j, k, grid, w)

    return qᶠᶠᶜ * wᶠᶠᶜ
end

function PVFluxX(mean_fields, sp)
    grid = mean_fields.u.grid
    loc = locationornothing((Center, Face, Face), mean_fields.u)
    return KernelFunctionOperation{loc...}(pv_flux_x_func, grid, sp, mean_fields)
end

function PVFluxY(mean_fields, sp)
    grid = mean_fields.u.grid
    loc = locationornothing((Face, Center, Face), mean_fields.u)
    return KernelFunctionOperation{loc...}(pv_flux_y_func, grid, sp, mean_fields)
end

function PVFluxZ(mean_fields, sp)
    grid = mean_fields.u.grid
    loc = locationornothing((Face, Face, Center), mean_fields.u)
    return KernelFunctionOperation{loc...}(pv_flux_z_func, grid, sp, mean_fields)
end
