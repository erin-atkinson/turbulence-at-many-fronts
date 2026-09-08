@inline function mld_func(i, j, _, grid, b, N²_min)
    k_below = grid.Nz + 1
    N²_below = ∂zᶜᶜᶠ(i, j, k_below, grid, b)
    
    for k in (grid.Nz + 1):-1:1
        k_below = k
        N²_below = ∂zᶜᶜᶠ(i, j, k_below, grid, b)
        N²_below >= N²_min && break
    end
    k_above = k_below + 1

    z_below, = node(i, j, k_below, grid, nothing, nothing, Face())
    z_above, = node(i, j, k_above, grid, nothing, nothing, Face())
    N²_above = ∂zᶜᶜᶠ(i, j, k_above, grid, b)

    k_below == 1 && return -z_below
    k_above == (grid.Nz + 1) && return -z_above

    z = z_below + (N²_min - N²_below) * (z_above - z_below) / (N²_above - N²_below)
    return -z
end

@doc raw"""
    MLD(b, N²_min)
Return a kernel function operation that calculates the mixed layer depth based on the first occurence of some buoyancy frequency.

See also [`MPEDensity`](@ref), [`COOLING`](@ref), [`MIXEDDensity`](@ref)

```math
N^2(x, y, -h_\text{ml}) = N^2_\text{min} \quad \text{where} \quad N^2 = \frac{\partial b}{\partial z}
```
"""
function MLD(b, N²_min)
    grid = b.grid
    loc = locationornothing((Center, Center, Nothing), b)
    return KernelFunctionOperation{loc...}(mld_func, grid, b, N²_min)
end