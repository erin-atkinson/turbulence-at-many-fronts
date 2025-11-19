@inline function N²_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    return ∂zᶜᶜᶠ(i, j, k, grid, dependency_fields.b_dfm)
end

N²_dependencies = (:b_dfm, )

@inline function M²_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    return ∂xᶠᶜᶜ(i, j, k, grid, dependency_fields.b_dfm)
end

M²_dependencies = (:b_dfm, )

@inline function Sx_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    return ℑxᶜᵃᵃ(i, j, k, grid, ∂zᶠᶜᶠ, dependency_fields.u_dfm)
end

Sx_dependencies = (:u_dfm, )

@inline function Sy_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    return ∂zᶜᶠᶠ(i, j, k, grid, dependency_fields.v_dfm)
end

Sx_dependencies = (:v_dfm, )

@inline function Rib_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    b_dfm = dependency_fields.b_dfm
    
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b_dfm)
    bz = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b_dfm)
    
    return sp.f^2 * bz / bx^2
end

Rib_dependencies = (:b_dfm, )

@inline function Ri_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    b_dfm = dependency_fields.b_dfm
    
    Sx = ℑxzᶜᵃᶜ(i, j, k, grid, ∂zᶠᶜᶠ, dependency_fields.u_dfm)
    Sy = ℑyzᵃᶜᶜ(i, j, k, grid, ∂zᶜᶠᶠ, dependency_fields.v_dfm)
    
    bz = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b_dfm)
    
    return bz / (Sx^2 + Sy^2)
end

Ri_dependencies = (:u_dfm, :v_dfm, :b_dfm)