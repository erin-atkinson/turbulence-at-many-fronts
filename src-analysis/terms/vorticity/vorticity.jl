@inline function ηx_func(i, j, k, grid, clock, sp, v, w)
    return ∂yᶜᶠᶠ(i, j, k, grid, w) - ∂zᶜᶠᶠ(i, j, k, grid, v)
end

@inline function ηy_func(i, j, k, grid, clock, sp, u, w)
    return ∂zᶠᶜᶠ(i, j, k, grid, u) - ∂xᶠᶜᶠ(i, j, k, grid, w)
end

@inline function ηz_func(i, j, k, grid, clock, sp, u, v)
    return ∂xᶠᶠᶜ(i, j, k, grid, v) - ∂yᶠᶠᶜ(i, j, k, grid, u) + sp.f
end
