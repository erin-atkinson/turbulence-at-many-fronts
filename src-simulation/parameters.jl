# parameters.jl

default_inputs = (;
    stop_time = 100e4, start_time = -40e4, save_time = 1e4, release_time = -1.0,
    f = 1e-4, L = 1e3,
    βx = 20, βh = 3,
    Nx = 1024, Nh = 768, Ny = 128, Nz = 64,
    βℓ = 1, βH = 0.1,
    βα = 0.1, βB = 0.0, βτ = 0.0, β₀ = 40, θτ = 0.0, βth = 5,
    wall_time_limit = 3600
)

@inline function create_simulation_parameters(input_parameters=(; ))
    ip = merge(default_inputs, input_parameters)

    fp = create_front_parameters(ip)

    # Domain size
    Lx = ip.βx * ip.L
    Lh = ip.βh * ip.L
    Ly = Lh * ip.Ny / ip.Nh
    Lz = 1.5fp.H

    # Strain
    α = ip.βα * ip.f

    # Cooling and mixing rate
    B = ip.βB * ip.L^2 * ip.f^3
    T_mix = (fp.H^2 / B)^(1/3)

    # Wind from turbulence scale?
    τ = ip.βτ * ip.L^2 * ip.f^2
    
    # Thermodynamics
    αV = 2.0678e-4 # K⁻¹
    cp = 4.1819e3 # J kg⁻¹ K⁻¹
    ρ = 1.027e3 # kg m⁻³
    g = 9.81 # m s⁻²
    
    Q_cool = (cp * ρ) * B / (αV * g)
    Q_stress = ρ * τ

    # Sponge layer damping rate
    σ = 0.5 * sqrt(fp.N₀²) / (2π)
    
    return merge(ip, fp, (; Lh, Lx, Ly, Lz, α, B, T_mix, τ, Q_cool, Q_stress, σ, αV, cp, ρ, g))
end

@inline function create_simulation_parameters(; input_parameters...)
    create_simulation_parameters(input_parameters)
end
