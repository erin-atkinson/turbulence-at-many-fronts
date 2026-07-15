@inline function buoyancy_density_func(i, j, k, grid, w, w_next, b)
    return @inbounds ℑzᵃᵃᶜ(i, j, k, grid, a_avg, w, w_next) * b[i, j, k]
end

"""
    BUOYANCYDensity(w, w_next, b)
Mean potential energy production
"""
function BUOYANCYDensity(w, w_next, b)
    grid = w.grid
    loc = locationornothing((Center, Center, Center), w)
    return KernelFunctionOperation{loc...}(buoyancy_density_func, grid, w, w_next, b)
end
