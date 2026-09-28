# The way conditional operators work is kinda coarse
# To future person: do this more intelligently
@inline function mixed_layer_mask_func(i, j, k, grid, h, ℓz)
    z = znode(k, grid, ℓz)
    Δz = zspacing(k, grid, ℓz)
    z_ml = @inbounds -h[i, j, k]
    
    z - Δz / 2 > z_ml && return one(eltype(grid))
    z + Δz / 2 < z_ml && return zero(eltype(grid))

    return 1 - (z_ml - (z - Δz / 2)) / Δz
end

@inline function mixed_layer_mask_func(i, j, k, grid, field, h, ℓz)
    val = @inbounds field[i, j, k]
    return mixed_layer_mask_func(i, j, k, grid, h, ℓz) * val
end

@inline function mixed_layer_average_func(i, j, k, grid, field, h, ℓz)
    return @inbounds mixed_layer_mask_func(i, j, k, grid, field, h, ℓz) / h[i, j, k]
end

@doc raw"""
    MixedLayerIntegral(field, h)
Return a kernel function operation that calculates an integral over the mixed layer

See also [`BuoyancyBasedMixedLayer`](@ref)
"""
function MixedLayerIntegral(field, h)
    grid = field.grid
    loc = location(field)
    return Integral(KernelFunctionOperation{loc...}(mixed_layer_mask_func, grid, field, h, loc[3]()); dims=3)
end

function MixedLayerAverage(field, h)
    grid = field.grid
    loc = location(field)
    return Integral(KernelFunctionOperation{loc...}(mixed_layer_average_func, grid, field, h, loc[3]()); dims=3)
end
