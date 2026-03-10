@inline function Rib_func(i, j, k, grid, clock, sp, b)
    bx = ℑxzᶜᵃᶠ(i, j, k, grid, ∂xᶠᶜᶜ, b)
    bz = ∂zᶜᶜᶠ(i, j, k, grid, b)
    
    return sp.f^2 * bz / bx^2
end

@inline function Ri_func(i, j, k, grid, clock, sp, u, v, b)
    Sx = ℑxᶜᵃᵃ(i, j, k, grid, ∂zᶠᶜᶠ, u)
    Sy = ℑyᵃᶜᵃ(i, j, k, grid, ∂zᶜᶠᶠ, v)
    bz = ∂zᶜᶜᶠ(i, j, k, grid, b)
    
    return bz / (Sx^2 + Sy^2)
end
