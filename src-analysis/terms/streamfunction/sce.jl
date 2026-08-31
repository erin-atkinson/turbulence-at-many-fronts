@inline function sce_density_func(i, j, k, grid, ψ)
    ∇²ψ = ∂xᶠᶜᶠ(i, j, k, grid, ∂xᶜᶜᶠ, ψ) + ∂zᶠᶜᶠ(i, j, k, grid, ∂zᶠᶜᶜ, ψ)
    
    return @inbounds -ψ[i, j, k] * ∇²ψ 
end

function SCEDensity(ψ)
    grid = ψ.grid
    loc = locationornothing((Face, Center, Face), ψ)
    return KernelFunctionOperation{loc...}(sce_density_func, grid, ψ)
end
