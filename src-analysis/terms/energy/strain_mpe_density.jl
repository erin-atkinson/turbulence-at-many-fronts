@inline function strain_mpe_density_func(i, j, k, grid, clock, b, h_ml, h_ml_next, sp)
    t = clock.time
    x, y, z = node(i, j, k, grid, Center(), Center(), Center())
    U = variable_strain_rate(t, sp) * velocity_profile(x, sp)

    h_ml_avg = a_avg(i, j, k, grid, h_ml, h_ml_next)
    
    return U * ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b) * (z + h_ml_avg)
end

"""
    STRAINMPEDensity(clock, b, h_ml, h_ml_next, sp)
"""
function STRAINMPEDensity(clock, b, h_ml, h_ml_next, sp)
    grid = b.grid
    loc = locationornothing((Center, Center, Center), b)
    return KernelFunctionOperation{loc...}(strain_mpe_density_func, grid, clock, b, h_ml, h_ml_next, sp)
end
