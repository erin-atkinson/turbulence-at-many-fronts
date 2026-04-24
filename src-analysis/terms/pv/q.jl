@inline function q_func(i, j, k, grid, sp, u, v, w, b)
    
    ηx = ℑxᶠᵃᵃ(i, j, k, grid, ηx_func, sp, v, w)
    ηy = ℑyᵃᶠᵃ(i, j, k, grid, ηy_func, sp, u, w)
    ηz = ℑzᵃᵃᶠ(i, j, k, grid, ηz_func, sp, u, v)

    bx = ∂xᶠᶜᶜ(i, j, k, grid, b)
    by = ∂yᶜᶠᶜ(i, j, k, grid, b)
    bz = ∂zᶜᶜᶠ(i, j, k, grid, b)
    
    return (bx * ηx
          + by * ηy
          + bz * ηz)
end

function PV(sp, u, v, w, b)
    grid = u.grid
    loc = locationornothing((Face, Face, Face), u)
    return KernelFunctionOperation{loc...}(q_func, grid, sp, u, v, w, b)
end