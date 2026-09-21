using Oceananigans.Operators
using Oceananigans.Grids: node
using CUDA: @allowscalar

# ---------------------------------------
# Quadratic damping mask halfway below thermocline
@inline function sponge_layer(x, y, z)
    H_sponge = sp.H + (sp.Lz - sp.H) / 2
    s = min((z+sp.Lz) / (sp.Lz - H_sponge), 1)
    return (1 - abs(s))^2
end
# ---------------------------------------

# ---------------------------------------
# Strain turns on slowly starting at t=0
@inline function variable_strain_rate(t, α, f, release_time)
    turnon = max(1-exp(-f * t / 10), 0)
    release_time <= 0 && return α * turnon
    
    turnoff = max(1-exp(-f * (t - release_time) / 10), 0)
    return α * (turnon - turnoff)
end

@inline function velocity_profile(x, Lh)
    return -2Lh * tanh(x / 2Lh)
end

# -d/dx velocity profile
@inline function strain_profile(x, Lh)
    return sech(x / 2Lh)^2
end

# Background velocity fields
const U = Field{Face, Nothing, Nothing}(grid)
const central_α = Field{Nothing, Nothing, Nothing}(grid)
const velocity_profile_array = Field{Face, Nothing, Nothing}(grid)
const strain_profile_array = Field{Center, Nothing, Nothing}(grid)

@allowscalar begin
    vp_func(x) = velocity_profile(x, sp.Lh)
    sp_func(x) = strain_profile(x, sp.Lh)
    
    velocity_profile_array.data.parent[:, 1, 1] .= vp_func.(grid.xᶠᵃᵃ.parent)
    strain_profile_array.data.parent[:, 1, 1] .= sp_func.(grid.xᶜᵃᵃ.parent)
end

# U is calculated directly to avoid issues with the boundary conditions
# As of this simulation, no open BCs for background fields like this
@inline function calculate_strain_callback(simulation, sp)
    t = simulation.model.clock.time
    α = variable_strain_rate(t, sp.α, sp.f, sp.release_time)
    
    set!(central_α, α)
    U.data.parent .= α .* velocity_profile_array.data.parent
    
    return nothing
end

@inline αf_func(i, j, k, grid, field) = @inbounds strain_profile_array[i, j, k] * central_α[i, j, k] * field[i, j, k]

@inline v_forcing_func(i, j, k, grid, clock, v) = -2αf_func(i, j, k, grid, v)
@inline w_forcing_func(i, j, k, grid, clock, w) = -αf_func(i, j, k, grid, w)
@inline b_forcing_func(i, j, k, grid, clock, b) = -αf_func(i, j, k, grid, b)
# ---------------------------------------

# ---------------------------------------
# Combination of forcings
u_forcing = (
    AdvectiveForcing(; u=U),
    Relaxation(; rate=sp.σ, mask=sponge_layer, target=0),
)
v_forcing = (
    AdvectiveForcing(; u=U),
    Relaxation(; rate=sp.σ, mask=sponge_layer, target=0),
    Forcing(v_forcing_func; discrete_form=true, field_dependencies=(:v, )),
)
w_forcing = (
    AdvectiveForcing(; u=U),
    Relaxation(; rate=sp.σ, mask=sponge_layer, target=0),
    Forcing(w_forcing_func; discrete_form=true, field_dependencies=(:w, )),
)
b_forcing = (
    AdvectiveForcing(; u=U),
    Forcing(b_forcing_func; discrete_form=true, field_dependencies=(:b, )),
)
# ---------------------------------------

forcing = (; u=u_forcing, v=v_forcing, w=w_forcing, b=b_forcing)
