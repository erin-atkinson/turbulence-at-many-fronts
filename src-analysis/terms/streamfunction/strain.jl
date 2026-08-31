@inline function strain_density_func(i, j, k, grid, clock, ψ, sp)
    x, = Oceananigans.Grids.node(i, j, k, grid, Face(), nothing, nothing)
    t = clock.time
    
    U = velocity_profile(x, sp) * variable_strain_rate(t, sp)
    ∇²ψ = ∂xᶠᶜᶠ(i, j, k, grid, ∂xᶜᶜᶠ, ψ) + ∂zᶠᶜᶠ(i, j, k, grid, ∂zᶠᶜᶜ, ψ)
    
    ψx = ℑxᶠᵃᵃ(i, j, k, grid, ∂xᶜᶜᶠ, ψ)

    return @inbounds -ψx * U * ∇²ψ
end

function STRAINDensity(clock, ψ, sp)
    grid = ψ.grid
    loc = locationornothing((Face, Center, Face), ψ)
    return KernelFunctionOperation{loc...}(strain_density_func, grid, clock, ψ, sp)
end
