#= simulation.jl
    Create a simulation of a front forced by strain flow and surface cooling

    Call with 
        julia -julia_opts -- simulation.jl ARGS...

    ARGS is a set of input arguments:
        [01]: Path to output folder
        [02]: Path to output folder to initialise from

        [03]: Simulation stop time in seconds
        [04]: Start time of simulation (negative for a cooling pre-initialisation)
        [05]: Strain release time
        [06]: Save interval

        [07]: Coriolis parameter in 1/seconds
        [08]: Deformation radius in metres

        [09]: Across-front domain size (relative to front)
        [10]: Central region domain size (relative to front)

        [11]: Number of across-front (x) grid points
        [12]: Number of central region (x) grid points
        [13]: Number of along-front (y) grid points
        [14]: Number of vertical (z) grid points

        [15]: Initial width of front relative to deformation radius
        [16]: Height of front relative to deformation radius

        [17]: Strain Rossby number
        [18]: Surface buoyancy flux 
        [19]: Surface wind stress
        [20]: Surface wind stress angle (0: along front pi/2: across front)
        [21]: Deep deformation radius
        [22]: Thermocline buoyancy step

        [23]: Maximum wall time
        
        [24]: Comment
=#
println("simulation.jl")
ENV["JULIA_SCRATCH_TRACK_ACCESS"] = 0
using Oceananigans

output_folder = ARGS[1]
initialised_folder = ARGS[2]

simulation_parameters = let

    # Simulation times
    stop_time = parse(Float64, ARGS[3])
    start_time = parse(Float64, ARGS[4])
    release_time = parse(Float64, ARGS[5])
    save_time = parse(Float64, ARGS[6])

    # Coriolis frequency
    f = parse(Float64, ARGS[7])

    # Deformation radius
    L = parse(Float64, ARGS[8])

    # Grid extent
    βx = parse(Float64, ARGS[9])
    βh = parse(Float64, ARGS[10])

    # Grid sizes
    Nx = parse(Int64, ARGS[11])
    Nh = parse(Int64, ARGS[12])
    Ny = parse(Int64, ARGS[13])
    Nz = parse(Int64, ARGS[14])

    # Front
    βℓ = parse(Float64, ARGS[15])
    βH = parse(Float64, ARGS[16])

    # Background
    βα = parse(Float64, ARGS[17])
    βB = parse(Float64, ARGS[18])
    βτ = parse(Float64, ARGS[19])
    θτ = parse(Float64, ARGS[20])
    β₀ = parse(Float64, ARGS[21])
    βth = parse(Float64, ARGS[22])

    wall_time_limit = parse(Float64, ARGS[23])
    
    (;
        stop_time, start_time, save_time, release_time,
        f, L,
        βx, βh,
        Nx, Nh, Ny, Nz,
        βℓ, βH,
        βα, βB, βτ, θτ, β₀, βth,
        wall_time_limit
    )
end

const comment = join(ARGS[24:end], " ")

for (k, v) in pairs(simulation_parameters)
    k_str = string(k)
    println("$k_str: $v")
end

include("create_simulation.jl")
