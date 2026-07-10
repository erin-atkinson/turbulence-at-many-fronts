#= simulation.jl
    Create a simulation of a front forced by strain flow and surface cooling

    Call with 
        julia -julia_opts -- simulation.jl ARGS...

    ARGS is a set of input arguments:
        [01]: Path to output folder

        [02]: Simulation stop time in seconds
        [03]: Start time of simulation (negative for a cooling pre-initialisation)
        [04]
        [05]: Save interval

        [06]: Coriolis parameter in 1/seconds
        [07]: Deformation radius in metres

        [08]: Across-front domain size (relative to front)
        [09]: Central region domain size (relative to front)

        [10]: Number of across-front (x) grid points
        [11]: Number of central region (x) grid points
        [12]: Number of along-front (y) grid points
        [13]: Number of vertical (z) grid points

        [14]: Initial width of front relative to deformation radius
        [15]: Height of front relative to deformation radius

        [16]: Strain Rossby number
        [17]: Surface buoyancy flux 
        [18]: Surface wind stress
        [19]: Surface wind stress angle (0: along front pi/2: across front)
        [20]: Deep deformation radius
        
        [21]: Comment
=#
ENV["JULIA_SCRATCH_TRACK_ACCESS"] = 0
using Oceananigans

output_folder = ARGS[1]

simulation_parameters = let

    # Simulation times
    stop_time = parse(Float64, ARGS[2])
    start_time = parse(Float64, ARGS[3])
    max_time = parse(Float64, ARGS[4])
    save_time = parse(Float64, ARGS[5])

    # Coriolis frequency
    f = parse(Float64, ARGS[6])

    # Deformation radius
    L = parse(Float64, ARGS[7])

    # Grid extent
    βx = parse(Float64, ARGS[8])
    βh = parse(Float64, ARGS[9])

    # Grid sizes
    Nx = parse(Int64, ARGS[10])
    Nh = parse(Int64, ARGS[11])
    Ny = parse(Int64, ARGS[12])
    Nz = parse(Int64, ARGS[13])

    # Front
    βℓ = parse(Float64, ARGS[14])
    βH = parse(Float64, ARGS[15])

    # Background
    βα = parse(Float64, ARGS[16])
    βB = parse(Float64, ARGS[17])
    βτ = parse(Float64, ARGS[18])
    θτ = parse(Float64, ARGS[19])
    β₀ = parse(Float64, ARGS[20])

    comment = join(ARGS[21:end], " ")
    
    (;
        stop_time, start_time, save_time, max_time,
        f, L,
        βx, βh,
        Nx, Nh, Ny, Nz,
        βℓ, βH,
        βα, βB, βτ, θτ, β₀,
        comment
    )
end

include("create_simulation.jl")
