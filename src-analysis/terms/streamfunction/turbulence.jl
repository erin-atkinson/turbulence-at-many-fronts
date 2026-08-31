@inline function turbulence_density_func(i, j, k, grid, ψ, turbulent_fluxes)
    uu = turbulent_fluxes.uu
    uw = turbulent_fluxes.uw
    wu = turbulent_fluxes.wu
    ww = turbulent_fluxes.ww
    
    𝔉uz = -∂zᶠᶜᶠ(i, j, k, grid, ∂xᶠᶜᶜ, uu) - ∂zᶠᶜᶠ(i, j, k, grid, ∂zᶠᶜᶜ, wu)
    𝔉wx = -∂xᶠᶜᶠ(i, j, k, grid, ∂xᶜᶜᶠ, uw) - ∂xᶠᶜᶠ(i, j, k, grid, ∂zᶜᶜᶠ, ww)
    
    𝔉ψ  = 𝔉uz - 𝔉wx
    
    return @inbounds ψ[i, j, k] * 𝔉ψ
end

function TURBULENCEDensity(ψ, turbulent_fluxes)
    grid = ψ.grid
    loc = locationornothing((Face, Center, Face), ψ)
    return KernelFunctionOperation{loc...}(turbulence_density_func, grid, ψ, turbulent_fluxes)
end
