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

"""
    MLD(b, N²_min)
Mixed layer depth via first instance of some buoyancy frequency
"""
function MLD(b, N²_min)
    grid = b.grid
    loc = locationornothing((Center, Center, Nothing), b)
    return KernelFunctionOperation{loc...}(mld_func, grid, b, N²_min)
end