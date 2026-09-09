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

@inline σff_avg(i, j, k, grid, loc, f, f_prev, sp) = f_avg_Gg(i, j, k, grid, f, f_prev, sponge_func, loc, f, sp)

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

@inline function b_flux_func(i, j, k, grid, clock, sp)
    t = clock.time
    return b_flux_func(t, sp)
end

# θ: angle relative to a down-front wind
# We only include wind in the central region
@inline function u_flux_func(x, t, sp) 
    turnon = 1 - exp(-sp.f*(t - sp.start_time) / 20)
    return -sp.τ * turnon * sin(sp.θτ) * exp(-x^2 / 4sp.L^2)
end

@inline function v_flux_func(x, t, sp) 
    turnon = 1 - exp(-sp.f*(t - sp.start_time) / 20)
    return -sp.τ * turnon * cos(sp.θτ) * exp(-x^2 / 4sp.L^2)
end

@inline function u_flux_func(i, j, k, grid, clock, sp)
    x, y, z = node(i, j, k, grid, Face(), Center(), Center())
    t = clock.time

    return u_flux_func(x, t, sp) 
end

@inline function v_flux_func(i, j, k, grid, clock, sp)
    x, y, z = node(i, j, k, grid, Center(), Face(), Center())
    t = clock.time

    return v_flux_func(x, t, sp) 
end

function UFlux(grid, clock, sp)
    return KernelFunctionOperation{Face, Center, Center}(u_flux_func, grid, clock, sp)
end

function VFlux(grid, clock, sp)
    return KernelFunctionOperation{Center, Face, Center}(v_flux_func, grid, clock, sp)
end

function BFlux(grid, clock, sp)
    return KernelFunctionOperation{Center, Center, Center}(b_flux_func, grid, clock, sp)
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

@inline function αff_avg(i, j, k, grid, loc, t, f, f_prev, sp)
    x, y, z = node(i, j, k, grid, loc...)

    α = variable_strain_rate(t, sp) * strain_profile(x, sp)
    f_avg  = a_avg(i, j, k, grid, f, f_prev)

    return @inbounds α * f[i, j, k] * f_avg
end
# ------------------------------------------------------------------------------
