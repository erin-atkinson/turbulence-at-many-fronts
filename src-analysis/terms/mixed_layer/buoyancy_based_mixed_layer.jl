@inline function buoyancy_based_mixed_layer_func(i, j, _, grid, b, δb)
    max_buoyancy = @inbounds b[i, j, 1]
    
    for k in 1:grid.Nz
        _b = @inbounds b[i, j, k]
        max_buoyancy = max(max_buoyancy, _b)
    end

    target_buoyancy = max_buoyancy - δb
    
    # We want the lowest occurence of this buoyancy, to ignore the inversion
    k_before = -1
    for k in 1:grid.Nz
        b_difference = @inbounds b[i, j, k] - target_buoyancy
        b_difference > 0 && break
        k_before = k
    end
    k_after = k_before + 1
    
    # Special escapes
    # Ensuring mixed layer is within domain just so calculations work later
    k_before == -1 && return -znode(1, grid, Face())
    k_before == grid.Nz && return -znode(grid.Nz, grid, Center())

    z_before = znode(k_before, grid, Center())
    z_after = znode(k_after, grid, Center())

    b_before = @inbounds b[i, j, k_before]
    b_after = @inbounds b[i, j, k_after]

    # Linear interpolation
    z = z_before + (target_buoyancy - b_before) * (z_after - z_before) / (b_after - b_before)
    
    # We return mixed layer _depths_
    return -z
end

@doc raw"""
    BuoyancyBasedMixedLayer(b, δb)
Return a kernel function operation that calculates the mixed layer depth based on the departure from the maximum buoyancy in a profile

See also [`MPEDensity`](@ref), [`COOLING`](@ref), [`MIXEDDensity`](@ref)
"""
function BuoyancyBasedMixedLayer(b, δb)
    grid = b.grid
    loc = locationornothing((Center, Center, Nothing), b)
    return KernelFunctionOperation{loc...}(buoyancy_based_mixed_layer_func, grid, b, δb)
end
