@inline function PE_func(i, j, k, grid, b)
    
    x, y, z = node(i, j, k, grid, Center(), Center(), Center())
    
    return @inbounds -b[i, j, k] * z
end

"""
    PE(b)
Potential energy
"""
function PE(b)
    grid = b.grid
    loc = locationornothing((Center, Center, Center), b)
    return KernelFunctionOperation{loc...}(PE_func, grid, b)
end
