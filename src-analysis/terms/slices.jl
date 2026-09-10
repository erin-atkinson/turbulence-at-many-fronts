using Oceananigans.Grids: node
using Oceananigans.Fields: interpolate

@inline function x_slice_func(i, j, k, grid, field, x, loc)
    (_, y, z) = node(i, j, k, grid, loc...)
    at_node = (x, y, z)

    return interpolate(at_node, field, loc, grid)
end

@inline function y_slice_func(i, j, k, grid, field, y, loc)
    (x, _, z) = node(i, j, k, grid, loc...)
    at_node = (x, y, z)
    
    return interpolate(at_node, field, loc, grid)
end

@inline function z_slice_func(i, j, k, grid, field, z, loc)
    (x, y, _) = node(i, j, k, grid, loc...)
    at_node = (x, y, z)
    
    return interpolate(at_node, field, loc, grid)
end

function XSlice(field, x)
    grid = field.grid
    T = eltype(grid)

    loc = location(field)
    instantiated_loc = instantiated_location(field)
    return KernelFunctionOperation{Nothing, loc[2], loc[3]}(x_slice_func, grid, field, T(x), instantiated_loc)
end

function YSlice(field, y)
    grid = field.grid
    T = eltype(grid)

    loc = location(field)
    instantiated_loc = instantiated_location(field)
    return KernelFunctionOperation{loc[1], Nothing, loc[3]}(y_slice_func, grid, field, T(y), instantiated_loc)
end

function ZSlice(field, z)
    grid = field.grid
    T = eltype(grid)

    loc = location(field)
    instantiated_loc = instantiated_location(field)
    return KernelFunctionOperation{loc[1], loc[2], Nothing}(z_slice_func, grid, field, T(z), instantiated_loc)
end