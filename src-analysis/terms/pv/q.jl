@inline function q_func(i, j, k, grid, clock, sp, u, v, w, b)
    
    ηx = ℑxᶠᵃᵃ(i, j, k, grid, ηx_func, clock, sp, v, w)
    ηy = ℑyᵃᶠᵃ(i, j, k, grid, ηy_func, clock, sp, u, w)
    ηz = ℑzᵃᵃᶠ(i, j, k, grid, ηz_func, clock, sp, u, v)

    bx = ∂xᶠᶜᶜ(i, j, k, grid, b)
    by = ∂yᶜᶠᶜ(i, j, k, grid, b)
    bz = ∂zᶜᶜᶠ(i, j, k, grid, b)
    
    return (bx * ηx
          + by * ηy
          + bz * ηz)
end

