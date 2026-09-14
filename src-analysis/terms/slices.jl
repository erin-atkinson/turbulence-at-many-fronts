# slices.jl
# Some functions for computing slices of fields along a direction
using Oceananigans.Utils: interpolator, _interpolate
using Oceananigans.Fields: fractional_x_index, fractional_y_index, fractional_z_index

@inline function x_slice_func(i, j, k, grid, field, x, loc)
    fractional_index =  fractional_x_index(x, loc, grid)
    ix = interpolator(fractional_index)
    data = @inbounds field[:, j, k]
    
    return _interpolate(data, ix)
end

@inline function y_slice_func(i, j, k, grid, field, y, loc)
    fractional_index =  fractional_y_index(y, loc, grid)
    iy = interpolator(fractional_index)
    data = @inbounds field[i, :, k]
    
    return _interpolate(data, iy)
end

@inline function z_slice_func(i, j, k, grid, field, z, loc)
    fractional_index =  fractional_z_index(z, loc, grid)
    iz = interpolator(fractional_index)
    data = @inbounds field[i, j, :]
    
    return _interpolate(data, iz)
end

@doc raw"""
    XSlice(field, x)
Return a KernelFunctionOperation that reduces an ND field to an N-1D slice at constant x
"""
function XSlice(field, x)
    grid = field.grid
    T = eltype(grid)

    loc = location(field)
    instantiated_loc = instantiated_location(field)
    return KernelFunctionOperation{Nothing, loc[2], loc[3]}(x_slice_func, grid, field, T(x), instantiated_loc)
end

@doc raw"""
    YSlice(field, y)
Return a KernelFunctionOperation that reduces an ND field to an N-1D slice at constant y
"""
function YSlice(field, y)
    grid = field.grid
    T = eltype(grid)

    loc = location(field)
    instantiated_loc = instantiated_location(field)
    return KernelFunctionOperation{loc[1], Nothing, loc[3]}(y_slice_func, grid, field, T(y), instantiated_loc)
end

@doc raw"""
    ZSlice(field, z)
Return a KernelFunctionOperation that reduces an ND field to an N-1D slice at constant z
"""
function ZSlice(field, z)
    grid = field.grid
    T = eltype(grid)

    loc = location(field)
    instantiated_loc = instantiated_location(field)
    return KernelFunctionOperation{loc[1], loc[2], Nothing}(z_slice_func, grid, field, T(z), instantiated_loc)
end
