using Oceananigans.Operators
using Oceananigans.Grids: node
using CUDA: @allowscalar

# ---------------------------------------
# Quadratic damping mask
@inline function sponge_layer(x, y, z)
    s = min((z+sp.Lz) / (sp.Lz-sp.H), 1)
    return (1 - abs(s))^2
end

# Damp b towards the bottom value
@inline function b_forcing_func(i, j, k, grid, clock, model_fields)
    (x, y, z, ) = node(i, j, k, grid, Center(), Center(), Center())
    (z_bottom, ) = node(i, j, 0, grid, Nothing(), Nothing(), Center())

    b = @inbounds model_fields.b[i, j, k]
    tb = @inbounds model_fields.b[i, j, 0] + sp.N₀² * (z - z_bottom)
    
    return sp.σ * min(tb - b, 0) * sponge_layer(x, y, z)
end
# ---------------------------------------

# ---------------------------------------
# Strain turns on slowly starting at t=0
@inline function variable_strain_rate(t, α, f)
    turnon = max(1-exp(-f * t / 15), 0)
    return α * turnon
end

@inline function velocity_profile(x, Lh)
    return -2Lh * tanh(x / 2Lh)
end

@inline function strain_profile(x, Lh)
    return sech(x / 2Lh)^2
end

# Background velocity fields
U = Field{Face, Nothing, Nothing}(grid)

# This is probably unnecessary
velocity_profile_array = Field{Face, Nothing, Nothing}(grid)
@allowscalar let vp(x) = velocity_profile(x, sp.Lh)
    velocity_profile_array.data.parent[:, 1, 1] .= vp.(grid.xᶠᵃᵃ.parent)
end

# U is calculated directly to avoid issues with boundary conditions
# Open never really worked...
@inline function calculate_U_callback(simulation, sp)
    t = simulation.model.clock.time
    α = variable_strain_rate(t, sp.α, sp.f)
    
    U.data.parent .= α .* velocity_profile_array.data.parent
    
    return nothing
end
# ---------------------------------------

# ---------------------------------------
# Background velocity forcing
@inline αf_func(x, y, z, t, f) = -variable_strain_rate(t, sp.α, sp.f) * strain_profile(x, sp.Lh) * f
@inline v_forcing_func(x, y, z, t, v) = 2αf_func(x, y, z, t, v)
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
    Forcing(v_forcing_func; field_dependencies=(:v, )),
)
w_forcing = (
    AdvectiveForcing(; u=U),
    Relaxation(; rate=sp.σ, mask=sponge_layer, target=0),
    Forcing(αf_func; field_dependencies=(:w, )),
)
b_forcing = (
    AdvectiveForcing(; u=U),
    Forcing(b_forcing_func; discrete_form=true),
    Forcing(αf_func; field_dependencies=(:b, )),
)
# ---------------------------------------

forcing = (; u=u_forcing, v=v_forcing, w=w_forcing, b=b_forcing)