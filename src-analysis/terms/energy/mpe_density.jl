@inline function mpe_density_func(i, j, k, grid, b, h_ml)
    
    x, y, z = node(i, j, k, grid, Center(), Center(), Center())
    
    return @inbounds -b[i, j, k] * (z + h_ml)
end

"""
    MPEDensity(b, h_ml)
Potential energy density
"""
function MPEDensity(b, h_ml)
    grid = b.grid
    loc = locationornothing((Center, Center, Center), b)
    return KernelFunctionOperation{loc...}(mpe_density_func, grid, b, h_ml)
end
