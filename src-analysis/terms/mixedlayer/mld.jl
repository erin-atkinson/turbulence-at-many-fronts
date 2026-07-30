using Oceananigans.Operators
using Oceananigans: location

# Integrate a field down to the mixed layer depth
@inline function ml_integrate_func(i, j, k, grid, ℓx, ℓy, ℓz, mld, field)
    zs = znodes(grid, ℓx, ℓy, ℓz)
    Δzs = zspacings(grid, ℓx, ℓy, ℓz)
    
    res = zero(eltype(field))
    h = @inbounds mld[i, j, k]
    
    k_above = findfirst(zs .> -h)
    k_above == nothing && return zero(eltype(field))
    k_below = max(k_above - 1, 1)

    z_above = zs[k_above]
    z_below = zs[k_below]

    δfield = @inbounds field[i, j, k_above] - field[i, j, k_below]
    δz = z_above - z_below
    δh = -h - z_below
    
    for _k in k_above:length(zs)
        res += @inbounds field[i, j, _k] * Δzs[_k]
    end
    
    k_above == k_below && return res
    
    if δh >= δz / 2
        res -= @inbounds field[i, j, k_above] * (δh - δz/2)
    else
        res += @inbounds field[i, j, k_below] * (δz/2 - δh)
    end
    
    return res
end

# Average a field within the mixed layer
@inline function ml_average_func(i, j, k, grid, ℓx, ℓy, ℓz, mld, field)
    return @inbounds ml_integrate_func(i, j, k, grid, ℓx, ℓy, ℓz, mld, field) / mld[i, j, k]
end


function ML_Average(field, mld)
    (ℓx, ℓy, ℓz) = location(field)
    grid = field.grid
    return KernelFunctionOperation{ℓx, ℓy, Nothing}(ml_average_func, grid, ℓx(), ℓy(), ℓz(), mld, field)
end

# Interpolate a field's value at the mixed layer depth
@inline function ml_interpolate_func(i, j, k, grid, ℓx, ℓy, ℓz, mld, field)
    zs = znodes(grid, ℓx, ℓy, ℓz)
    h = @inbounds mld[i, j, k]
    
    k_above = findfirst(zs .> -h)
    
    k_above == nothing && return @inbounds field[i, j, length(zs)]
    k_above == 1 && return @inbounds field[i, j, 1]
    
    k_below = k_above - 1

    z_above = zs[k_above]
    z_below = zs[k_below]

    δfield = @inbounds field[i, j, k_above] - field[i, j, k_below]
    δz = z_above - z_below
    δh = -h - z_below
    
    return @inbounds field[i, j, k_below] + δfield / δz * δh
end

function ML_Interpolate(field, mld)
    (ℓx, ℓy, ℓz) = location(field)
    grid = field.grid
    
    return KernelFunctionOperation{ℓx, ℓy, Nothing}(ml_interpolate_func, grid, ℓx(), ℓy(), ℓz(), mld, field)
end

@inline function mld_func(i, j, _, grid, b, ε)
    # Mixed layer depth is depth at which the b differs
    # from the surface value by ε
    zs = znodes(grid, nothing, nothing, Center())
    Δzs = zspacings(grid, nothing, nothing, Center())
    
    # average from _k to the surface
    b_monotonic = zeros(eltype(b), grid.Nz)
    b_monotonic[1] = @inbounds b[i, j, 1]
    
    for k in 2:grid.Nz
        Δb = @inbounds max(b[i, j, k] - b[i, j, k-1], zero(eltype(b)))
        b_monotonic[k] = b_monotonic[k-1] + Δb
    end

    b_surface = b_monotonic[grid.Nz]
    b_mld = b_surface - ε

    k_above = findfirst(b_monotonic .> b_mld)
    
    k_above == nothing && return -zs[grid.Nz]
    k_above == 1 && return -zs[1]
    
    k_below = max(k_above - 1, 1)

    z_above = zs[k_above]
    z_below = zs[k_below]

    δb = b_monotonic[k_above] - b_monotonic[k_below]
    δz = z_above - z_below

    h = -(z_below + δz / δb * (b_mld - b_monotonic[k_below]))
    return h
end

@inline function constant_mld_func(i, j, k, grid, b, ε)
    res, = node(i, j, 1, grid, nothing, nothing, Center())
    res = -res
    return res
    for i in axes(b, 1), j in axes(b, 2)
        res = min(mld_func(i, j, k, grid, b, ε), res)
    end
    return res
end

function MLD(b, ε)
    (ℓx, ℓy, ℓz) = location(b)
    grid = b.grid
    KernelFunctionOperation{ℓx, ℓy, Nothing}(mld_func, grid, b, ε)
end

function ConstantMLD(b, ε)
    (ℓx, ℓy, ℓz) = location(b)
    grid = b.grid
    KernelFunctionOperation{Nothing, Nothing, Nothing}(constant_mld_func, grid, b, ε)
end

# The tolerance required to get an mld of H at t=0
function ε_estimate(initial_b, sp)
    zs = znodes(initial_b)
    b_profile = mean(interior(initial_b, :, :, :); dims=(1, 2))[1, 1, :]
    k_above = findfirst(zs .> -sp.H)
    k_below = k_above - 1
    
    return mean(b_profile[k_below:end]) - b_profile[k_below]
end

function growth_rate_func(i, j, k, grid, clock, mld, mld_prev)
    mld_avg = @inbounds (mld[i, j, k] + mld_prev[i, j, k]) / 2
    mld_diff = @inbounds (mld[i, j, k] - mld_prev[i, j, k])
    return mld_diff / clock.last_Δt / mld_avg
end

function GrowthRate(clock, mld, mld_prev)
    (ℓx, ℓy, ℓz) = location(mld)
    grid = mld.grid
    KernelFunctionOperation{ℓx, ℓy, ℓz}(growth_rate_func, grid, clock, mld, mld_prev)
end