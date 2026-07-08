using Oceananigans.Fields: ZeroField

# Strain turns on slowly starting at t=0
@inline function variable_strain_rate(t, α, f, max_time)
    turnon = max(1-exp(-f * t / 15), 0)
    max_time <= 0 && return α * turnon
    
    turnoff = max(1-exp(-f * (t - max_time) / 15), 0)
    return α * (turnon - turnoff)
end

@inline function velocity_profile(x, Lh)
    return -2Lh * tanh(x / 2Lh)
end

@inline function strain_profile(x, Lh)
    return sech(x / 2Lh)^2
end

@inline variable_strain_rate(t) = variable_strain_rate(t, sp.α, sp.f, sp.max_time)
@inline velocity_profile(x) = velocity_profile(x, sp.Lh)
@inline strain_profile(x) = strain_profile(x, sp.Lh)

# Background velocity fields
U = Field{Face, Nothing, Nothing}(grid)
V = ZeroField()
W = ZeroField()

# This is probably unnecessary
velocity_profile_array = Field{Face, Nothing, Nothing}(grid)
velocity_profile_array.data.parent[:, 1, 1] .= velocity_profile.(grid.xᶠᵃᵃ.parent)

# U and V are calculated directly to avoid issues with boundary conditions
# Open does NOT work
@inline function compute_background!(U, V, W, clock)
    t = clock.time
    α = variable_strain_rate(t)
    U.data.parent .= α .* velocity_profile_array.data.parent
    return nothing
end
