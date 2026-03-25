using Oceananigans.Grids: xnode, ynode, znode
using Oceananigans: location

struct AbstractKernel end

"""
    coarse_grain(i, j, k, grid, kernel, loc, field)
Kernel function operation that returns a coarse-grained field
"""
for ξ in (:x, :y, :z)
    coarse_grain_ξ = Symbol(:coarse_grain, ξ)
    @eval begin
        $coarse_grain_ξ(i, j, k, grid, kernel, loc, field) = @inbounds field[i, j, k]
        $coarse_grain_ξ(i, j, k, grid, kernel, loc, f::Function, args...) = f(i, j, k, grid, args...)
    end
end

@inline f(i, j, k, grid, weight, field) = @inbounds field[i, j, k] * weight
@inline f(i, j, k, grid, weight, func::Function, args...) = func(i, j, k, grid, args...) * weight

const XPeriodicGrid = AbstractGrid{<:Any, Periodic}
const YPeriodicGrid = AbstractGrid{<:Any, <:Any, Periodic}
const ZPeriodicGrid = AbstractGrid{<:Any, <:Any, <:Any, Periodic}

@inline function ingrid(i, j, k, grid, ℓx, ℓy, ℓz)
    Nx = grid.Nx + ℓx isa Face
    Ny = grid.Ny + ℓy isa Face
    Nz = grid.Nz + ℓz isa Face

    x = 1 <= i <= Nx || grid isa XPeriodicGrid
    y = 1 <= j <= Ny || grid isa YPeriodicGrid
    z = 1 <= k <= Nz || grid isa ZPeriodicGrid

    return x && y && z
end

function coarse_grain_x(i, j, k, grid, loc, kernel::AbstractKernel, args...)
    hw = half_width_x(i, j, k, grid, loc, kernel)
    x = xnode(i, j, k, grid, loc, nothing, nothing)

    res = 0.0
    weights = 0.0

    for di in -hw:hw
        !ingrid(i + di, j, k, grid, loc, nothing, nothing) && continue

        dx = xnode(i + di, j, k, grid, loc, nothing, nothing) - x
        weight = kernel_func(kernel, dx)
        res += f(i + di, j, k, grid, weight, args...)
        weights += weight
    end

    return res / weight
end

function coarse_grain_y(i, j, k, grid, loc, kernel::AbstractKernel, field)
    hw = half_width_y(i, j, k, grid, loc, kernel)
    y = ynode(i, j, k, grid, nothing, loc, nothing)

    res = 0.0
    weights = 0.0

    for dj in -hw:hw
        !ingrid(i, j + dj, k, grid, nothing, loc, nothing) && continue

        dy = ynode(i, j + dj, k, grid, nothing, loc, nothing) - y
        weight = kernel_func(kernel, dy)
        res += f(i, j + dj, k, grid, weight, args...)
        weights += weight
    end

    return res / weight
end

function coarse_grain_z(i, j, k, grid, loc, kernel::AbstractKernel, field)
    hw = half_width_z(i, j, k, grid, loc, kernel)
    z = znode(i, j, k, grid, nothing, nothing, loc)

    res = 0.0
    weights = 0.0

    for dk in -hw:hw
        !ingrid(i, j, k + dk, grid, nothing, nothing, loc) && continue

        dz = znode(i, j, k + dk, grid, nothing, nothing, loc) - z
        weight = kernel_func(kernel, dz)
        res += f(i, j, k + dk, grid, weight, args...)
        weights += weight
    end

    return res / weight
end

function coarse_grain(i, j, k, grid, ℓx, ℓy, ℓz, kernels, field)
    coarse_grain_x(i, j, k, grid, ℓx, kernels[1],
        coarse_grain_y, ℓy, kernels[2], 
        coarse_grain_z, ℓz, kernels[3], 
        field
    )
end

function Coarse(field, k1, k2, k3)
    (ℓx, ℓy, ℓz) = location(field)
    kernels = (k1, k2, k3)
    grid = field.grid
    
    return KernelFunctionOperation{ℓx, ℓy, ℓz}(coarse_grain, grid, ℓx(), ℓy(), ℓz(), kernels, field)
end

"""
    struct Gaussian{S}
Gaussian coarse-graining kernel
"""
struct Gaussian{S} <: AbstractKernel
    σ :: S
end

@inline kernel_func(kernel::Gaussian, dx) = exp(-dx^2 / 2kernel.σ^2)

@inline function half_width_x(i, j, k, grid, loc, kernel::Gaussian)
    return Integer(1 + (5 * kernel.σ) / minimum_xspacing(grid))
end

@inline function half_width_y(i, j, k, grid, loc, kernel::Gaussian)
    return Integer(1 + (5 * kernel.σ) / minimum_yspacing(grid))
end

@inline function half_width_z(i, j, k, grid, loc, kernel::Gaussian)
    return Integer(1 + (5 * kernel.σ) / minimum_zspacing(grid))
end
