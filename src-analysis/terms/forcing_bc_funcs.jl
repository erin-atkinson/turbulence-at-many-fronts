# Sponge layer -----------------------------------------------------------------
# Quadratic damping mask
@inline function sponge_layer_func(z, sp)
    s = min((z+sp.Lz) / (sp.Lz-sp.H), 1)
    return sp.σ * (1 - abs(s))^2
end

@inline function sponge_func(i, j, k, grid, loc, field, sp)
    z = znode(i, j, k, grid, loc...)
    σ = sponge_layer_func(z, sp)
    return @inbounds -σ * field[i, j, k]
end

@inline function sponge_func(i, j, k, grid, loc, f, sp, args...)
    z = znode(i, j, k, grid, loc...)
    σ = sponge_layer_func(z, sp)
    return -σ * f(i, j, k, grid, args...)
end

# Damp b towards the bottom value
@inline function b_forcing_func(i, j, k, grid, b, sp)
    (x, y, z, ) = node(i, j, k, grid, Center(), Center(), Center())
    (z_bottom, ) = node(i, j, 1, grid, Nothing(), Nothing(), Center())

    b = @inbounds b[i, j, k]
    tb = @inbounds b[i, j, 1] + sp.N₀² * (z - z_bottom)
    
    return sp.σ * min(tb - b, 0) * sponge_layer_func(z, sp)
end

function SpongeLayer(field, sp)
    (ℓx, ℓy, ℓz) = location(field)
    grid = field.grid
    
    return KernelFunctionOperation{ℓx, ℓy, ℓz}(sponge_func, grid, (ℓx(), ℓy(), ℓz()), field, sp)
end
# ------------------------------------------------------------------------------

# Surface ----------------------------------------------------------------------
# Cooling turns on slowly
@inline function b_flux_func(t, sp) 
    turnon = 1 - exp(-sp.f*(t - sp.start_time) / 20)
    return sp.B * turnon
end

# θ: angle relative to a down-front wind
# We only include wind in the central region
@inline function u_flux_func(x, y, t, sp) 
    turnon = 1 - exp(-sp.f*(t - sp.start_time) / 20)
    return -sp.τ * turnon * sin(sp.θτ) * exp(-x^2 / 4sp.L^2)
end

@inline function v_flux_func(x, y, t, sp) 
    turnon = 1 - exp(-sp.f*(t - sp.start_time) / 20)
    return -sp.τ * turnon * cos(sp.θτ) * exp(-x^2 / 4sp.L^2)
end
# ------------------------------------------------------------------------------

# Background flow --------------------------------------------------------------
@inline function variable_strain_rate(t, sp)
    turnon = max(1-exp(-sp.f * t / 15), 0)
    sp.max_time <= 0 && return sp.α * turnon
    
    turnoff = max(1-exp(-sp.f * (t - sp.max_time) / 15), 0)
    return sp.α * (turnon - turnoff)
end

@inline function velocity_profile(x, sp)
    return -2sp.Lh * tanh(x / 2sp.Lh)
end

@inline function strain_profile(x, sp)
    return sech(x / 2sp.Lh)^2
end
# ------------------------------------------------------------------------------
