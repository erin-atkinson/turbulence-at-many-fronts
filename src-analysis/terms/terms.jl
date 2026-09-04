using Oceananigans.Operators
using Oceananigans.Grids: node
using Oceananigans: location

@inline along_front_mean(a) = Field(Average(a; dims=2))
@inline afm(a) = along_front_mean(a)

@inline fg(i, j, k, grid, f, g) = @inbounds f[i, j, k] * g[i, j, k]
@inline fGg(i, j, k, grid, f, G, args...) = @inbounds f[i, j, k] * G(i, j, k, grid, args...)
@inline FfGg(i, j, k, grid, F, f, G, args...) = @inbounds F(i, j, k, grid, f) * G(i, j, k, grid, args...)

@inline f_avg(i, j, k, grid, f, f_next) = @inbounds (f[i, j, k] + f_next[i, j, k]) / 2

@inline f_avg_Gg(i, j, k, grid, f, f_next, G, args...) = f_avg(i, j, k, grid, f, f_next) * G(i, j, k, grid, args...)

locationornothing(loc, u) = map(loc, location(u)) do ℓ, ℓu
    ℓu isa Type{Nothing} ? ℓu : ℓ
end

include("CoarseGraining.jl")
include("constants.jl")

# Vorticity and gradients
include("gradients/vorticity.jl")
include("gradients/richardson.jl")
