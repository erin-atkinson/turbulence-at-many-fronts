@inline function bflux_density_func(i, j, k, grid, wb)
    return ℑzᵃᵃᶜ(i, j, k, grid, wb)
end

"""
    BFLUXDensity(wb)
Buoyancy flux
"""
function BFLUXDensity(wb)
    grid = wb.grid
    loc = locationornothing((Center, Center, Center), wb)
    return KernelFunctionOperation{loc...}(bflux_density_func, grid, wb)
end
