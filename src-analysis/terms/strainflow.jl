using Oceananigans.Fields: ZeroField

# Strain turns on slowly starting at t=0
@inline variable_strain_rate(t) = variable_strain_rate(t, sp)
@inline velocity_profile(x) = velocity_profile(x, sp)
@inline strain_profile(x) = strain_profile(x, sp)

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
