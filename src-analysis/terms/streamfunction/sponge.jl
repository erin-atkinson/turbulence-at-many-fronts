
@inline function sponge_density_func(i, j, k, grid, ψ, sp)
    Suz = -∂zᶠᶜᶠ(i, j, k, grid, sponge_func, (Center(), nothing, Face()), ∂zᶠᶜᶜ, sp, ψ)
    Swx = ∂xᶠᶜᶠ(i, j, k, grid, sponge_func, (Face(), nothing, Center()), ∂xᶜᶜᶠ, sp, ψ)
    
    return Suz - Swx
end

function SPONGEDensity(ψ, sp)
    grid = ψ.grid
    loc = locationornothing((Face, Center, Face), ψ)
    return KernelFunctionOperation{loc...}(sponge_density_func, grid, ψ, sp)
end
