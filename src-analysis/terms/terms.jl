using Oceananigans.Operators
using Oceananigans.Grids: node

@inline down_front_mean(a) = Field(Average(a; dims=2))
@inline dfm(a) = down_front_mean(a)

# Calculate a down front-mean of a field
# Grid spacing is uniform in down-front direction
@inline function discrete_down_front_mean(i, j, k, grid, f, args...)
    j_indices = (grid.Hy+1):(grid.Hy+grid.Ny)
    total = mapreduce(+, j_indices) do _j
        f(i, _j, k, grid, args...)
    end
    total / grid.Ny
end

@inline function discrete_down_front_mean(i, j, k, grid, field)
    j_indices = (grid.Hy+1):(grid.Hy+grid.Ny)
    total = mapreduce(+, j_indices) do _j
        @inbounds field[i, _j, k]
    end
    total / grid.Ny
end

@inline ddfm(args...) = discrete_down_front_mean(args...)

@inline fg(i, j, k, grid, f, g) = @inbounds f[i, j, k] * g[i, j, k]
@inline fGg(i, j, k, grid, f, G, args...) = @inbounds f[i, j, k] * G(i, j, k, grid, args...)
@inline FfGg(i, j, k, grid, F, f, G, args...) = @inbounds F(i, j, k, grid, f) * G(i, j, k, grid, args...)

@inline a_avg(i, j, k, grid, a, a_next) = @inbounds (a[i, j, k] + a_next[i, j, k]) / 2

@inline f_avg_Gg(i, j, k, grid, f, f_next, G, args...) = a_avg(i, j, k, grid, f, f_next) * G(i, j, k, grid, args...)
@inline f′_avg_Gg(i, j, k, grid, f, f_next, f_dfm, f_next_dfm, G, args...) = (a_avg(i, j, k, grid, f, f_next) - a_avg(i, j, k, grid, f_dfm, f_next_dfm)) * G(i, j, k, grid, args...)

@inline df′dt(i, j, k, grid, f, f_next, f_dfm, f_next_dfm, Δt) = (f′(i, j, k, grid, f_next, f_next_dfm) - f′(i, j, k, grid, f, f_dfm)) / Δt

@inline f′(i, j, k, grid, f, f_dfm) = @inbounds f[i, j, k] - f_dfm[i, j, k]
@inline f′g′(i, j, k, grid, f, f_dfm, g, g_dfm) = f′(i, j, k, grid, f, f_dfm) * f′(i, j, k, grid, g, g_dfm)
@inline f′Gg′(i, j, k, grid, f, f_dfm, G, g, g_dfm) = f′(i, j, k, grid, f, f_dfm) * G(i, j, k, grid, f′, g, g_dfm)

#=
# Faster to assume that the kernel is separable and normalized
@inline function coarse_grain_x(i, j, k, grid, weights, f::F, args...) where {F <: Function}
    res = 0
    for d in axes(weights, 1)
        w = @inbounds weights[d]
        ind =  min(max(i+d - size(weights, 1) ÷ 2, 1), grid.Nx)
        res += f(ind, j, k, grid, args...) * w
    end

    return res
end

@inline function coarse_grain_y(i, j, k, grid, weights, f::F, args...) where {F <: Function}
    res = 0
    for d in axes(weights, 1)
        w = @inbounds weights[d]
        ind =  min(max(j+d - size(weights, 1) ÷ 2, 1), grid.Ny)
        res += f(i, ind, k, grid, args...) * w
    end

    return res
end

@inline function coarse_grain_z(i, j, k, grid, weights, f::F, args...) where {F <: Function}
    res = 0
    for d in axes(weights, 1)
        w = @inbounds weights[d]
        ind =  min(max(k+d - size(weights, 1) ÷ 2, 1), grid.Nz)
        res += f(i, j, ind, grid, args...) * w
    end

    return res
end

@inline identityc(i, j, k, grid, field) = @inbounds field[i, j, k]
@inline coarse_grain_x(i, j, k, grid, weights, field) = coarse_grain_x(i, j, k, grid, weights, identityc, field)
@inline coarse_grain_y(i, j, k, grid, weights, field) = coarse_grain_y(i, j, k, grid, weights, identityc, field)
@inline coarse_grain_z(i, j, k, grid, weights, field) = coarse_grain_z(i, j, k, grid, weights, identityc, field)

@inline function coarse_grain_variable_x(i, j, k, grid, Lx, kernel_size, σ, args...)
    weights = gaussian_weights_x(i, j, k, grid, Lx, kernel_size, σ)
    return coarse_grain_x(i, j, k, grid, weights, args...)
end

# 1D array of normalized weights
@inline function gaussian_weights(Δx, σ)
    weights = map(Δx) do x
        exp(-x^2 / 2σ^2)
    end
    return weights ./ sum(weights)
end

# 1D array of normalized weights
@inline function gaussian_weights_x(i, j, k, grid, Lx, kernel_size, σ)
    x, = node(i, j, k, grid, Lx, Nothing(), Nothing())

    Δx = map(-kernel_size:kernel_size) do di
        i2 = min(max(i+di, 1), grid.Nx)
        x1, = node(i2, j, k, grid, Lx, Nothing(), Nothing())
        x1 - x
    end

    return gaussian_weights(Δx, σ)
end
=#

include("strainflow.jl")
include("CoarseGraining.jl")
include("constants.jl")