# parameters.jl

default_inputs = (;
    run_time = 8e5, start_time = -4e5, save_time = 1e3,
    f = 1e-4, L = 1e3,
    βx = 4, βh = 1,
    Nx = 1024, Nh = 768, Ny = 128, Nz = 256,
    βb = 1, βℓ = 1, βH = 0.1,
    α = 1e-5, Q = 100, β₀ = 8,
    comment = ""
)

@inline function create_simulation_parameters(input_parameters=(; ))
    ip = merge(default_inputs, input_parameters)

    fp = create_front_parameters(ip)

    # Domain size
    Lx = ip.βx * ip.L
    Lh = ip.βh * ip.L
    Ly = Lh * ip.Ny / ip.Nh
    Lz = 1.5fp.H

    # Physical constants needed to construct surface BC
    αV = 2.0678e-4 # K⁻¹
    cₚ = 4.1819e3 # J kg⁻¹ K⁻¹
    ρ = 1.027e3 # kg m⁻³
    g = 9.81 # m s⁻²

    # Surface buoyancy flux
    B = αV * g * ip.Q / (cₚ * ρ)

    # Sponge layer damping rate
    σ = 0.5 * sqrt(fp.N₀²) / (2π)
    
    return merge(ip, fp, (; Lh, Lx, Ly, Lz, B, σ))
end

@inline function create_simulation_parameters(; input_parameters...)
    create_simulation_parameters(input_parameters)
end
