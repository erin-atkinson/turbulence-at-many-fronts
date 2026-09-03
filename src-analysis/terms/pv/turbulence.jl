@inline function turbulent_pv_flux_x_func(i, j, k, grid, sp, mean_fields, turbulent_forcing)
    # x-flux of PV: cff
    𝔉vᶜᶠᶠ = ℑzᵃᵃᶠ(i, j, k, grid, turbulent_forcing.v)
    𝔉wᶜᶠᶠ = ℑyᵃᶠᵃ(i, j, k, grid, turbulent_forcing.w)
    𝔇ᶜᶠᶠ = ℑyzᵃᶠᶠ(i, j, k, grid, turbulent_forcing.b)

    ηx  = ηx_func(i, j, k, grid, sp, mean_fields.v, mean_fields.w)
    
    byᶜᶠᶠ = ℑzᵃᵃᶠ(i, j, k, grid, ∂yᶜᶠᶜ, mean_fields.b)
    bzᶜᶠᶠ = ℑyᵃᶠᵃ(i, j, k, grid, ∂zᶜᶜᶠ, mean_fields.b)
    
    return -(𝔉vᶜᶠᶠ * bzᶜᶠᶠ - 𝔉wᶜᶠᶠ * byᶜᶠᶠ) - ηx * 𝔇ᶜᶠᶠ
end

@inline function turbulent_pv_flux_y_func(i, j, k, grid, sp, mean_fields, turbulent_forcing)
    # y-flux of PV: fcf
    𝔉uᶠᶜᶠ = ℑzᵃᵃᶠ(i, j, k, grid, turbulent_forcing.u)
    𝔉wᶠᶜᶠ = ℑxᶠᵃᵃ(i, j, k, grid, turbulent_forcing.w)
    𝔇ᶠᶜᶠ = ℑxzᶠᵃᶠ(i, j, k, grid, turbulent_forcing.b)

    ηy  = ηy_func(i, j, k, grid, sp, mean_fields.u, mean_fields.w)
    
    bxᶠᶜᶠ = ℑzᵃᵃᶠ(i, j, k, grid, ∂xᶠᶜᶜ, mean_fields.b)
    bzᶠᶜᶠ = ℑxᶠᵃᵃ(i, j, k, grid, ∂zᶜᶜᶠ, mean_fields.b)
    
    return -(𝔉wᶠᶜᶠ * bxᶠᶜᶠ - 𝔉uᶠᶜᶠ * bzᶠᶜᶠ) - ηy * 𝔇ᶠᶜᶠ
end

@inline function turbulent_pv_flux_z_func(i, j, k, grid, sp, mean_fields, turbulent_forcing)
    # z-flux of PV: ffc
    𝔉uᶠᶠᶜ = ℑyᵃᶠᵃ(i, j, k, grid, turbulent_forcing.u)
    𝔉vᶠᶠᶜ = ℑxᶠᵃᵃ(i, j, k, grid, turbulent_forcing.v)
    𝔇ᶠᶠᶜ = ℑxyᶠᶠᵃ(i, j, k, grid, turbulent_forcing.b)

    ηz  = ηz_func(i, j, k, grid, sp, mean_fields.u, mean_fields.v)
    
    bxᶠᶠᶜ = ℑyᵃᶠᵃ(i, j, k, grid, ∂xᶠᶜᶜ, mean_fields.b)
    byᶠᶠᶜ = ℑxᶠᵃᵃ(i, j, k, grid, ∂yᶜᶠᶜ, mean_fields.b)
    
    return -(𝔉uᶠᶠᶜ * byᶠᶠᶜ - 𝔉vᶠᶠᶜ * bxᶠᶠᶜ) - ηz * 𝔇ᶠᶠᶜ
end

function TurbulentPVFluxX(mean_fields, turbulent_forcing, sp)
    grid = mean_fields.u.grid
    loc = locationornothing((Center, Face, Face), mean_fields.u)
    return KernelFunctionOperation{loc...}(turbulent_pv_flux_x_func, grid, sp, mean_fields, turbulent_forcing)
end

function TurbulentPVFluxY(mean_fields, turbulent_forcing, sp)
    grid = mean_fields.u.grid
    loc = locationornothing((Face, Center, Face), mean_fields.u)
    return KernelFunctionOperation{loc...}(turbulent_pv_flux_y_func, grid, sp, mean_fields, turbulent_forcing)
end

function TurbulentPVFluxZ(mean_fields, turbulent_forcing, sp)
    grid = mean_fields.u.grid
    loc = locationornothing((Face, Face, Center), mean_fields.u)
    return KernelFunctionOperation{loc...}(turbulent_pv_flux_z_func, grid, sp, mean_fields, turbulent_forcing)
end
