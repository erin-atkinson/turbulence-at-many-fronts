using JLD2
using Oceananigans.TimeSteppers: Clock

function check_completion(simname)
    foldername = "../scratch/turbulence-at-many-fronts/$simname"
    
    !isdir(foldername) && return "$simname: Uninitialised (no folder)"
    filenames = readdir(foldername)
    !mapreduce(a->startswith(a, "checkpoint"), |, filenames; init=false) && return "$simname: Uninitialised (no checkpoint)"
    
    checkpointnames = filter(a->startswith(a, "checkpoint"), filenames)
    prev_time = mapreduce(max, checkpointnames; init=-Inf) do checkpoint_file
        str = "simulation/model/clock"
        checkpoint_path = joinpath(foldername, checkpoint_file)
        
        jldopen(file->file[str].time, checkpoint_path)
    end
    prev_time < 0 && return "$simname: Uninitialised (no checkpoint at t=0)"
    
    !mapreduce(a->startswith(a, "OUTPUT"), |, filenames) && return "$simname: Initialised only (no output)"
    
    iterations = jldopen(joinpath(foldername, "OUTPUT.jld2")) do file
        keys(file["timeseries/t"])
    end
    ts = jldopen(joinpath(foldername, "OUTPUT.jld2")) do file
        [file["timeseries/t/$iteration"] for iteration in keys(file["timeseries/t"])]
    end
    sp = jldopen(joinpath(foldername, "OUTPUT.jld2")) do file
        file["metadata/parameters"]
    end
    "$simname: Run until ft = $(sp.f * ts[end]), αt = $(sp.α * ts[end]), $(iterations[end]) $checkpointnames"
end

function make_filename(sp, ext=nothing, pre=""; βH=sp.βH, βα=sp.βα, βB=sp.βB, βτ=sp.βτ, θτ=sp.θτ)
    strs = map([βH, βα, βB, βτ, θτ]) do β
        replace(string(β), "."=>"_")
    end
    ext = ext == nothing ? "" : ".$ext"
    return joinpath(pre, join(strs, "-") * ext)
end

function θτ_from_str(ip)
    ip.βτ == 0 && return "C"
    ip.θτ == 0 && return "N"
    ip.θτ == π/2 && return "E"
end

function make_preamble(jobname, T)
    """
    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --gpus-per-node=1
    #SBATCH --time=$T
    #SBATCH --job-name=$jobname
    #SBATCH --output=../scratch/logs/$jobname.txt
    
    module load julia/1.12.5
    
    # Launch from scratch
    export JULIA_DEPOT_PATH=\$SCRATCH/julia-trig
    # export JULIA_CUDA_SOFT_MEMORY_LIMIT=10%
    cd ~/turbulence-at-many-fronts
    """
end

function make_body(ip, filename)
    """
    
    output_folder=\$SCRATCH/turbulence-at-many-fronts/$filename

    stop_time=$(ip.stop_time)
    start_time=$(ip.start_time)
    max_time=$(ip.max_time)
    save_time=$(ip.save_time)
    
    f=$(ip.f)
    L=$(ip.L)
    
    betax=$(ip.βx)
    betah=$(ip.βh)
    
    Nx=$(ip.Nx)
    Nh=$(ip.Nh)
    Ny=$(ip.Ny)
    Nz=$(ip.Nz)
    
    betal=$(ip.βℓ)
    betaH=$(ip.βH)
    
    betaalpha=$(ip.βα)
    betaB=$(ip.βB)
    betatau=$(ip.βτ)
    thetatau=$(ip.θτ)
    beta0=$(ip.β₀)
    
    comment="$(ip.comment)"

    julia -t 8 -- src-simulation/simulation.jl \$output_folder \$stop_time \$start_time \$max_time \$save_time \$f \$L \$betax \$betah \$Nx \$Nh \$Ny \$Nz \$betal \$betaH \$betaalpha \$betaB \$betatau \$thetatau \$beta0 \$comment
    """
end

function make_copy(filename::String)
    """
    
    mkdir \$output_folder/../$filename
    cp \$output_folder/checkpoint* \$output_folder/../$filename
    """
end

function make_copy(filenames)
    mapreduce(make_copy, *, filenames; init="")
end


function make_script(jobname, T, ip, filename, destfilenames=[])
    preamble = make_preamble(jobname, T)
    body = make_body(ip, filename)
    return preamble * body * make_copy(destfilenames)
end

function save_script(jobname, T, ip, filename; loc="", copy_to=[])
    scriptpath = joinpath(loc, jobname * ".sh")
    println(check_completion(filename))
    write(scriptpath, make_script(jobname, T, ip, filename, copy_to))
    return nothing
end

include("ensemble.jl")

# All the simulations run until wall-time is 3 hours, then checkpoint.
# Run this file again to check for completion, then submit others
T = "4:00:00"

println()
println("TEST: VARYING CENTRAL REGION SIZE")
path = "jobs-simulation/size-test"
mkpath(path)

for (ip, filename) in zip(size_test.ips, size_test.filenames)
    save_script(filename, T, ip, filename; loc=path)
end

println()
println("TEST: VARYING RESOLUTION")
path = "jobs-simulation/resolution-test"
mkpath(path)

for (ip, filename) in zip(resolution_test.ips, resolution_test.filenames)
    save_script(filename, T, ip, filename; loc=path)
end

println()
println("TEST: VARYING OUTER REGION RESOLUION")
path = "jobs-simulation/outer-test"
mkpath(path)

for (ip, filename) in zip(outer_test.ips, outer_test.filenames)
    save_script(filename, T, ip, filename; loc=path)
end

println()
println("VARYING COOLING AND DEPTH: INITIALISATION")
path = "jobs-simulation/cooling-depth-init"
mkpath(path)

for (ip, filename, copy_to) in zip(cooling_depth_init.ips, cooling_depth_init.filenames, cooling_depth_init.destfilenamess)
    save_script(filename, T, ip, filename; loc=path, copy_to)
end

println()
println("VARYING COOLING, DEPTH AND STRAIN")
path = "jobs-simulation/cooling-depth"
mkpath(path)

println("01")
for (ip, filename) in zip(cooling_depth_01.ips, cooling_depth_01.filenames)
    save_script(filename, T, ip, filename; loc=path)
end

println("02")
for (ip, filename) in zip(cooling_depth_02.ips, cooling_depth_02.filenames)
    save_script(filename, T, ip, filename; loc=path)
end

println("005")
for (ip, filename) in zip(cooling_depth_005.ips, cooling_depth_005.filenames)
    save_script(filename, T, ip, filename; loc=path)
end

#=
println("Cooling only")
for (ip, filename) in zip(cooling_only.ips, cooling_only.filenames)
    save_script(filename, T, ip, filename; loc=path)
end

println("Depth only")
for (ip, filename) in zip(depth_only.ips, depth_only.filenames)
    save_script(filename, T, ip, filename; loc=path)
end
=#
# Varying along-front winds

# Varying across-front winds
