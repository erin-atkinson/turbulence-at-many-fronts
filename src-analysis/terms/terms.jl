using Oceananigans.Operators
using Oceananigans.Grids: node
using Oceananigans: location, instantiated_location

@inline along_front_mean(a) = Field(Average(a; dims=2))
@inline afm(a) = along_front_mean(a)

@inline fg(i, j, k, grid, f, g) = @inbounds f[i, j, k] * g[i, j, k]
@inline fGg(i, j, k, grid, f, G, args...) = @inbounds f[i, j, k] * G(i, j, k, grid, args...)
@inline FfGg(i, j, k, grid, F, f, G, args...) = @inbounds F(i, j, k, grid, f) * G(i, j, k, grid, args...)

@inline f_avg(i, j, k, grid, f, f_prev) = @inbounds (f[i, j, k] + f_prev[i, j, k]) / 2
@inline f_avg_Gg(i, j, k, grid, f, f_prev, G, args...) = f_avg(i, j, k, grid, f, f_prev) * G(i, j, k, grid, args...)

locationornothing(loc, u) = map(loc, location(u)) do ℓ, ℓu
    ℓu isa Type{Nothing} ? ℓu : ℓ
end

include("CoarseGraining.jl")
include("constants.jl")
include("slices.jl")
include("forcing_bc_funcs.jl")

# Vorticity and gradients
include("gradients/vorticity.jl")
include("gradients/richardson.jl")

# Helpers for advection terms
include("advection/advection.jl")
include("advection/diffusion.jl")
include("advection/operators.jl")

# Mean potential energy
include("energy/mld.jl")
include("energy/mpe_density.jl")
include("energy/bflux_density.jl")
include("energy/mixed_density.jl")
include("energy/cooling.jl")
include("energy/strain_mpe_density.jl")

# Mean kinetic
include("energy/mke_density.jl")
include("energy/dsp_density.jl")
include("energy/stress.jl")
include("energy/buoyancy_density.jl")
include("energy/sponge_mke_density.jl")
include("energy/lsp_density.jl")
include("energy/vsp_density.jl")
include("energy/strain_mke_density.jl")

# Potential vorticity
include("potential_vorticity/potential_vorticity.jl")
include("potential_vorticity/flux.jl")
include("potential_vorticity/background.jl")
include("potential_vorticity/turbulence.jl")

# Streamfunction energy
include("streamfunction/streamfunction.jl")
include("streamfunction/sce_density.jl")
include("streamfunction/ageostrophic_density.jl")
include("streamfunction/turbulence_sce_density.jl")
include("streamfunction/sponge_sce_density.jl")
include("streamfunction/strain_sce_density.jl")
