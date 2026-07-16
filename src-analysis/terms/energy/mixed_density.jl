@inline function mixed_density_func(i, j, k, grid, clock, b, b_next, h_ml, h_ml_next)
    Δt = clock.last_Δt

    b_avg = a_avg(i, j, k, grid, b, b_next)
    ∂h_ml∂t = @inbounds (h_ml_next[i, j, k] - h_ml[i, j, k]) / Δt

    return -b_avg * ∂h_ml∂t
end

"""
    MIXEDDensity(clock, b, b_next, h_ml, h_ml_next)
"""
function MIXEDDensity(clock, b, b_next, h_ml, h_ml_next)
    grid = b.grid
    loc = locationornothing((Center, Center, Center), b)
    return KernelFunctionOperation{loc...}(mixed_density_func, grid, clock, b, b_next, h_ml, h_ml_next)
end
