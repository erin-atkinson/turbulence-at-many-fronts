# base_state.jl
# Functions describing the initial state of the simulations

using Oceananigans: fill_halo_regions!

@inline f(s) = (1 + tanh(s)) / 2 # 0 to 1
@inline f′(s) = (sech(s)^2) / 2
const f′′_max = maximum(s->-tanh(s) * sech(s)^2, range(-10, 10, 1000))

@inline g(s) = log(1 + exp(s))
@inline g′(s) = 1 / (1 + exp(-s))
@inline g′′(s) = exp(-s) / (1 + exp(-s))^2

# Background stratification
@inline b∞(z, sp) = -sp.λ * sp.H * sp.N₀² * g(-(z + sp.H) / (sp.λ * sp.H))

# Buoyancy
@inline function front_buoyancy(x, z, sp)
    x₁ = x / sp.ℓ + sp.a * (z + sp.H / 2) / sp.H
    z₁ = (z + sp.H) / (sp.λ * sp.H)
    
    return sp.Δb * f(x₁) * g′(z₁) + b∞(z, sp)
end

# Stratification
@inline function front_N²(x, z, sp)
    x₁ = x / sp.ℓ + sp.a * (z + sp.H / 2) / sp.H
    z₁ = (z + sp.H) / (sp.λ * sp.H)
    
    return (
          sp.a * (sp.Δb / sp.H) * f′(x₁) * g′(z₁)
        + (sp.Δb / (sp.λ * sp.H)) * f(x₁) * g′′(z₁)
        + sp.N₀² * g′(-z₁)
    )
end

# Horizontal buoyancy gradient
@inline function front_M²(x, z, sp)
    x₁ = x / sp.ℓ + sp.a * (z + sp.H / 2) / sp.H
    z₁ = (z + sp.H) / (sp.λ * sp.H)
    
    return (sp.Δb / sp.ℓ) * f′(x₁) * g′(z₁)
end

# Thermal wind shear
@inline front_S(x, z, sp) = front_M²(x, z, sp) / sp.f

@inline function approximate_front_velocity(x, z, sp)
    x₁ = x / sp.ℓ + sp.a * (z + sp.H / 2) / sp.H
    x₂ = x / sp.ℓ - sp.a / 2
    z₁ = (z + sp.H) / (sp.λ * sp.H)
    
    return (sp.H * sp.Δb / (sp.a * sp.ℓ * sp.f)) * (f(x₁) - f(x₂)) * g′(z₁)
end

@inline function approximate_surface_velocity(x, sp)
    x₁ = x / sp.ℓ + sp.a / 2
    
    return (sp.H * sp.Δb / (sp.ℓ * sp.f)) * f′(x₁)
end

@inline front_Ri(x, z, sp) = front_N²(x, z, sp) / front_shear(x, z, sp)^2

@inline function maximum_front_velocity(ip)
    return ip.L * ip.f / ip.βℓ
end

@inline function maximum_Ro(ip)
    # Maximum Rossby number occurs at the surface
    return 2f′′_max * ip.f / ip.βℓ^2
end

@inline function create_front_parameters(ip)
    λ = 0.03

    # Buoyancy change of Ri = 1 front
    Δb = 2ip.L * ip.f^2 / ip.βH

    # Size relative to deformation radius
    ℓ = ip.βℓ * ip.L

    # MLD
    H = ip.βH * ip.L

    # Deep water buoyancy frequency from deformation radius
    N₀² = (ip.f * ip.β₀ * ip.L / H)^2

    Ro = maximum_Ro(ip)
    V = maximum_front_velocity(ip)
    a = 1/ip.βℓ^2
    
    return (; λ, Δb, ℓ, H, N₀², Ro, V, a)
end
create_front_parameters(; ip...) = create_front_parameters(ip)

@inline function front_initial_conditions(grid::RectilinearGrid, sp)
    # Use Oceananigans fields to setup the initial thermal wind properly
    
    b = Field{Center, Center, Center}(grid)
    set!(b, (x, y, z)->front_buoyancy(x, z, sp))
    fill_halo_regions!(b)

    # Compute the thermal wind shear
    S_op = @at (Center, Face, Face) ∂x(b) / sp.f
    S = compute!(Field(S_op))

    # Integrate
    V_op = CumulativeIntegral(S; dims=3)
    v = compute!(Field(V_op))
    fill_halo_regions!(v)
    
    # Random secondary circulation
    u(x, y, z) = 1e-14 * randn()
    w(x, y, z) = 1e-14 * randn()
    
    return (; u, v, w, b)
end
