@inline function ageostrophic_density_func(i, j, k, grid, ψ, v, b, sp)
    bx = ℑzᵃᵃᶠ(i, j, k, grid, ∂xᶠᶜᶜ, b)
    vz = ℑxᶠᵃᵃ(i, j, k, grid, ∂zᶜᶜᶠ, v)

    return @inbounds ψ[i, j, k] * (sp.f * vz - bx)
end

function AGEOSTROPHICDensity(ψ, v, b, sp)
    grid = ψ.grid
    loc = locationornothing((Face, Center, Face), ψ)
    return KernelFunctionOperation{loc...}(ageostrophic_density_func, grid, ψ, v, b, sp)
end
