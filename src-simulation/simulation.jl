#= simulation.jl
    Create a simulation of a front forced by strain flow and surface cooling

    Call with 
        julia -julia_opts -- simulation.jl ARGS...

    ARGS is a set of input arguments:
        [01]: Path to output folder

        [02]: Simulation run time in seconds
        [03]: Start time of simulation (negative for a cooling pre-initialisation)
        [04]: Save interval

        [05]: Coriolis parameter in 1/seconds
        [06]: Deformation radius in metres

        [07]: Across-front domain size (relative to front)
        [08]: Central region domain size (relative to front)

        [09]: Number of across-front (x) grid points
        [10]: Number of central region (x) grid points
        [11]: Number of along-front (y) grid points
        [12]: Number of vertical (z) grid points

        [13]: Initial width of front relative to deformation radius
        [14]: Height of front relative to deformation radius

        [15]: Strain Rossby number
        [16]: Surface buoyancy flux 
        [17]: Surface wind stress
        [18]: Surface wind stress angle (0: along front pi/2: across front)
        [19]: Deep deformation radius
        
        [20]: Comment
=#
ENV["JULIA_SCRATCH_TRACK_ACCESS"] = 0
using Oceananigans

output_folder = ARGS[1]

simulation_parameters = let

    # Simulation times
    run_time = parse(Float64, ARGS[2])
    start_time = parse(Float64, ARGS[3])
    save_time = parse(Float64, ARGS[4])

    # Coriolis frequency
    f = parse(Float64, ARGS[5])

    # Deformation radius
    L = parse(Float64, ARGS[6])

    # Grid extent
    βx = parse(Float64, ARGS[7])
    βh = parse(Float64, ARGS[8])

    # Grid sizes
    Nx = parse(Int64, ARGS[9])
    Nh = parse(Int64, ARGS[10])
    Ny = parse(Int64, ARGS[11])
    Nz = parse(Int64, ARGS[12])

    # Front
    βℓ = parse(Float64, ARGS[13])
    βH = parse(Float64, ARGS[14])

    # Background
    βα = parse(Float64, ARGS[15])
    βB = parse(Float64, ARGS[16])
    βτ = parse(Float64, ARGS[17])
    θτ = parse(Float64, ARGS[18])
    β₀ = parse(Float64, ARGS[19])

    comment = join(ARGS[20:end], " ")
    
    (;
        run_time, start_time, save_time,
        f, L,
        βx, βh,
        Nx, Nh, Ny, Nz,
        βℓ, βH,
        βα, βB, βτ, θτ, β₀,
        comment
    )
end

include("create_simulation.jl")
