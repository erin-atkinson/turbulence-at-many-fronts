using Oceananigans.Grids: xnode, ynode, znode, AbstractGrid
using Oceananigans: location

abstract type AbstractKernel end

for ξ in (:x, :y, :z)
    coarse_grain_ξ = Symbol(:coarse_grain_, ξ)
    @eval begin
        $coarse_grain_ξ(i, j, k, grid, loc, kernel, field) = @inbounds field[i, j, k]
        $coarse_grain_ξ(i, j, k, grid, loc, kernel, f::Function, args...) = f(i, j, k, grid, args...)
    end
end

@inline f(i, j, k, grid, weight, field) = @inbounds field[i, j, k] * weight
@inline f(i, j, k, grid, weight, func::Function, args...) = func(i, j, k, grid, args...) * weight

const XPeriodicGrid = AbstractGrid{<:Any, Periodic}
const YPeriodicGrid = AbstractGrid{<:Any, <:Any, Periodic}
const ZPeriodicGrid = AbstractGrid{<:Any, <:Any, <:Any, Periodic}

@inline function ingrid(i, j, k, grid, ℓx, ℓy, ℓz)
    Nx = grid.Nx + (ℓx isa Face)
    Ny = grid.Ny + (ℓy isa Face)
    Nz = grid.Nz + (ℓz isa Face)

    x = 1 <= i <= Nx || (grid isa XPeriodicGrid)
    y = 1 <= j <= Ny || (grid isa YPeriodicGrid)
    z = 1 <= k <= Nz || (grid isa ZPeriodicGrid)

    return x && y && z
end

function coarse_grain(i, j, k, grid, ℓx, ℓy, ℓz, kernel, field)
    σx = kernel.σx
    σy = kernel.σy
    σz = kernel.σz
    wx = kernel.wx
    wy = kernel.wy
    wz = kernel.wz
    
    (x, y, z) = node(i, j, k, grid, ℓx, ℓy, ℓz)
    
    res = 0.0
    weights = 0.0

    for di in -wx:wx
        _i, _j, _k = i + di, j, k
        
        !ingrid(_i, _j, _k, grid, ℓx, ℓy, ℓz) && continue
        (dx, dy, dz) = node(_i, _j, _k, grid, ℓx, ℓy, ℓz) .- (x, y, z)
        
        weight = kernel_func(kernel, dx, dy, dz)
        
        res += @inbounds field[_i, _j, _k] * weight
        weights += weight
    end

    for dj in -wy:wy
        _i, _j, _k = i, j + dj, k
        
        !ingrid(_i, _j, _k, grid, ℓx, ℓy, ℓz) && continue
        (dx, dy, dz) = node(_i, _j, _k, grid, ℓx, ℓy, ℓz) .- (x, y, z)
        
        weight = kernel_func(kernel, dx, dy, dz)
        
        res += @inbounds field[_i, _j, _k] * weight
        weights += weight
    end

    for dk in -wz:wz
        _i, _j, _k = i, j, k + dk
        
        !ingrid(_i, _j, _k, grid, ℓx, ℓy, ℓz) && continue
        (dx, dy, dz) = node(_i, _j, _k, grid, ℓx, ℓy, ℓz) .- (x, y, z)
        
        weight = kernel_func(kernel, dx, dy, dz)
        
        res += @inbounds field[_i, _j, _k] * weight
        weights += weight
    end
    
    weights == 0.0 && return 0.0
    return res / weights
end

regularize_location(loc) = loc
regularize_location(loc::Nothing) = Center()

function Coarse(field, kernel)
    (ℓx, ℓy, ℓz) = location(field)
    
    grid = field.grid
    locs = (regularize_location(a()) for a in (ℓx, ℓy, ℓz))
    return KernelFunctionOperation{ℓx, ℓy, ℓz}(coarse_grain, grid, locs..., kernel, field)
end

"""
    struct Gaussian{S}
Gaussian coarse-graining kernel
"""
struct Gaussian{G, S, I} <: AbstractKernel
    grid::G
    σx::S
    σy::S
    σz::S
    wx::I
    wy::I
    wz::I
end

function Gaussian(grid, σx, σy, σz)
    wx = Integer((5 * σx) ÷ minimum_xspacing(grid))
    wy = Integer((5 * σy) ÷ minimum_yspacing(grid))
    wz = Integer((5 * σz) ÷ minimum_zspacing(grid))

    return Gaussian(grid, σx, σy, σz, wx, wy, wz)
end

@inline kernel_func(kernel::Gaussian, dx, dy, dz) = exp(-dx^2 / 2kernel.σx^2) * exp(-dy^2 / 2kernel.σy^2) * exp(-dz^2 / 2kernel.σz^2)
