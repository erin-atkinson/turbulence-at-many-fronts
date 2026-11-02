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
    cd ~/turbulence-at-many-fronts
    """
end

function make_body(ip, filename, init_folder="")
    """
    
    output_folder=\$SCRATCH/turbulence-at-many-fronts/$filename
    init_folder=\$SCRATCH/turbulence-at-many-fronts/$init_folder
    
    stop_time=$(ip.stop_time)
    start_time=$(ip.start_time)
    release_time=$(ip.release_time)
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
    thetatau=$(float(ip.θτ))
    beta0=$(ip.β₀)
    betath=$(ip.βth)

    wall_time_limit=$(ip.wall_time_limit)
    
    comment="$(ip.comment)"

    julia -t 8 -- src-simulation/simulation.jl \$output_folder \$init_folder \$stop_time \$start_time \$release_time \$save_time \$f \$L \$betax \$betah \$Nx \$Nh \$Ny \$Nz \$betal \$betaH \$betaalpha \$betaB \$betatau \$thetatau \$beta0 \$betath \$wall_time_limit \$comment
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


function make_script(jobname, T, ip, filename; copy_to=[], init_folder="")
    preamble = make_preamble(jobname, T)
    body = make_body(ip, filename, init_folder)
    return preamble * body * make_copy(copy_to)
end

function save_script(jobname, T, ip, filename; loc="", copy_to=[], init_folder="")
    scriptpath = joinpath(loc, jobname * ".sh")
    println(check_completion(filename))
    write(scriptpath, make_script(jobname, T, ip, filename; copy_to, init_folder))
    return nothing
end

function make_mps_script(mpsjobname, jobnames, T; loc)
    jobs = ["\nbash " * joinpath(loc, jobname * ".sh") * " &\nsleep 10\n" for jobname in jobnames]
    
    """
    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --gpus-per-node=1
    #SBATCH --time=$T
    #SBATCH --job-name=$mpsjobname
    #SBATCH --output=../scratch/logs/$mpsjobname.txt
    
    export JULIA_CUDA_MEMORY_POOL=split
    export JULIA_CUDA_SOFT_MEMORY_LIMIT=10%
    rm -rf /dev/shm/$mpsjobname
    
    mkdir /dev/shm/$mpsjobname
    mkdir /dev/shm/$mpsjobname/nvidia-mps
    mkdir /dev/shm/$mpsjobname/nvidia-log
    
    export CUDA_MPS_PIPE_DIRECTORY=/dev/shm/$mpsjobname/nvidia-mps
    export CUDA_MPS_LOG_DIRECTORY=/dev/shm/$mpsjobname/nvidia-log
    export CUDA_VISIBLE_DEVICES=0
    
    nvidia-cuda-mps-control -d
    
    cd ~/turbulence-at-many-fronts
    
    module load julia/1.12.5
    export JULIA_DEPOT_PATH=\$SCRATCH/julia-trig
    $(join(jobs))
    pkill -c -f nvidia-cuda-mps-server
    wait
    
    echo quit | nvidia-cuda-mps-control
    """
end

function save_mps_script(mpsjobname, jobnames, T; loc="")
    scriptpath = joinpath(loc, mpsjobname * ".sh")
    write(scriptpath, make_mps_script(mpsjobname, jobnames, T; loc))
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

path = "jobs-simulation/cooling-only-init"
mkpath(path)

for (ip, filename, copy_to) in zip(cooling_only_init.ips, cooling_only_init.filenames, cooling_only_init.destfilenamess)
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

println("03")
for (ip, filename) in zip(cooling_depth_03.ips, cooling_depth_03.filenames)
    save_script(filename, T, ip, filename; loc=path)
end

println("05")
for (ip, filename) in zip(cooling_depth_05.ips, cooling_depth_05.filenames)
    save_script(filename, T, ip, filename; loc=path)
end

println("07")
for (ip, filename) in zip(cooling_depth_07.ips, cooling_depth_07.filenames)
    save_script(filename, T, ip, filename; loc=path)
end

println()
println("COOLING ONLY")
path = "jobs-simulation/cooling-only"
mkpath(path)

for (ip, filename) in zip(cooling_only.ips, cooling_only.filenames)
    save_script(filename, T, ip, filename; loc=path)
end

println()
println("VARYING WIND: INITIALISATION")
path = "jobs-simulation/wind-init"
mkpath(path)

for (ip, filename, copy_to) in zip(wind_init.ips, wind_init.filenames, wind_init.destfilenamess)
    save_script(filename, T, ip, filename; loc=path, copy_to)
end
save_mps_script("all-wind", wind_init.filenames[1:3], T; loc=path)

println()
println("VARYING WIND")
path = "jobs-simulation/wind"
mkpath(path)

for (ip, filename) in zip(wind.ips, wind.filenames)
    save_script(filename, T, ip, filename; loc=path)
end

println()
println("HIGH RESOLUTION: INITIALISATION")
path = "jobs-simulation/high-resolution-init"
mkpath(path)

T = "12:00:00"

for (ip, filename, copy_to, init_folder) in zip(highresolution_init.ips, highresolution_init.filenames, highresolution_init.destfilenamess, hd_init_folders)
    save_script(filename, T, ip, filename; loc=path, copy_to, init_folder)
end

println()
println("HIGH RESOLUTION")
path = "jobs-simulation/high-resolution"
mkpath(path)

for (ip, filename) in zip(highresolution.ips, highresolution.filenames)
    save_script(filename, T, ip, filename; loc=path)
end
