using JLD2

function check_completion(simname)
    foldername = "../scratch/turbulence-at-many-fronts/$simname"
    
    !isdir(foldername) && return "$simname: Uninitialised (no folder)"
    filenames = readdir(foldername)
    !mapreduce(a->startswith(a, "checkpoint"), |, filenames) && return "$simname: Uninitialised (no checkpoint)"
    !mapreduce(a->startswith(a, "OUTPUT"), |, filenames) && return "$simname: Initialised only (no output)"
    checkpointnames = filter(a->startswith(a, "checkpoint"), filenames)
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

    run_time=$(ip.run_time)
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

    julia -t 8 -- src-simulation/simulation.jl \$output_folder \$run_time \$start_time \$max_time \$save_time \$f \$L \$betax \$betah \$Nx \$Nh \$Ny \$Nz \$betal \$betaH \$betaalpha \$betaB \$betatau \$thetatau \$beta0 \$comment
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

default_ip = (;
    run_time = 8e5, start_time = -2e5, save_time = 3e3, max_time = -1.0,
    f = 1e-4, L = 1e3,
    βx = 20, βh = 3,
    Nx = 1024, Nh = 768, Ny = 128, Nz = 64,
    βℓ = 6, βH = 0.1,
    βα = 0.1, βB = 0.03, βτ = 0.00, θτ = 0.0, β₀ = 6,
    comment = ""
)
# All the simulations run until wall-time is 3 hours, then checkpoint.
# Run this file again to check for completion, then submit others
# Let's leave some wiggle room
T = "4:00:00"

println()
println("TEST: VARYING CENTRAL REGION SIZE")
path = "jobs-simulation/size-test"
mkpath(path)
ips = [
    (; default_ip..., βh=2, save_time=1e4, comment="Central region size test 2"),
    (; default_ip..., βh=3, save_time=1e4, comment="Central region size test 3"),
    (; default_ip..., βh=4, save_time=1e4, comment="Central region size test 4"),
    (; default_ip..., βh=5, save_time=1e4, comment="Central region size test 5"),
]
filenames = [
    "size-test-h2",
    "size-test-h3",
    "size-test-h4",
    "size-test-h5",
]
for (ip, filename) in zip(ips, filenames)
    save_script(filename, T, ip, filename; loc=path)
end

println()
println("TEST: VARYING RESOLUTION")
path = "jobs-simulation/resolution-test"
mkpath(path)
ips = [
    (; default_ip..., save_time=1e4, βh = 4, Nx=1400, Nh=1024, Nz=96, comment="Resolution test x1400 z96"),
    (; default_ip..., save_time=1e4, βh = 4, Nx=1400, Nh=1024, Nz=64, comment="Resolution test x1400 z64"),
    (; default_ip..., save_time=1e4, βh = 4, Nx=1024, Nh=768, Nz=96, comment="Resolution test x1024 z96"),
    (; default_ip..., save_time=1e4, βh = 4, Nx=1024, Nh=768, Nz=32, comment="Resolution test x1024 z32"),
    (; default_ip..., save_time=1e4, βh = 4, Nx=800, Nh=600, Nz=64, comment="Resolution test x800 z64"),
    (; default_ip..., save_time=1e4, βh = 4, Nx=800, Nh=600, Nz=32, comment="Resolution test x800 z32"),
]
filenames = [
    "resolution-test-x1400-z96",
    "resolution-test-x1400-z64",
    "resolution-test-x1024-z96",
    "resolution-test-x1024-z32",
    "resolution-test-x800-z64",
    "resolution-test-x800-z32",
]

for (ip, filename) in zip(ips, filenames)
    save_script(filename, T, ip, filename; loc=path)
end

println()
println("TEST: VARYING OUTER REGION RESOLUION")
path = "jobs-simulation/outer-test"
mkpath(path)
ips = [
    (; default_ip..., save_time=1e4, βh = 3, Nx=800, Nh=700, Nz=64, comment="Outer test x800"),
    (; default_ip..., save_time=1e4, βh = 3, Nx=900, Nh=700, Nz=64, comment="Outer test x900"),
    (; default_ip..., save_time=1e4, βh = 3, Nx=1000, Nh=700, Nz=64, comment="Outer test x1000"),
    (; default_ip..., save_time=1e4, βh = 3, Nx=1100, Nh=700, Nz=64, comment="Outer test x1100"),
]
filenames = [
    "outer-test-x800",
    "outer-test-x900",
    "outer-test-x1000",
    "outer-test-x1100",
]
for (ip, filename) in zip(ips, filenames)
    save_script(filename, T, ip, filename; loc=path)
end

println()
println("VARYING COOLING AND DEPTH: INITIALISATION")
default_ip = (;
    run_time = 8e5, start_time = -2e5, save_time = 1e4, max_time = -1.0,
    f = 1e-4, L = 1e3,
    βx = 20, βh = 3,
    Nx = 1400, Nh = 1024, Ny = 128, Nz = 64,
    βℓ = 6, βH = 0.1,
    βα = 0.1, βB = 0.03, βτ = 0.00, θτ = 0.0, β₀ = 6,
    comment = "Cooling, depth initialisation"
)

path = "jobs-simulation/cooling-depth-init"
mkpath(path)
ips = [
    (; default_ip..., Nz=39, βH=0.06, βB=0.01, run_time=4.3e5, start_time=-4.3e5),
    (; default_ip..., Nz=39, βH=0.06, βB=0.03, run_time=3e5, start_time=-3e5),
    (; default_ip..., Nz=39, βH=0.06, βB=0.05, run_time=2.5e5, start_time=-2.5e5),
    (; default_ip..., Nz=39, βH=0.06, βB=0.1, run_time=2e5, start_time=-2e5),
    (; default_ip..., Nz=39, βH=0.06, βB=0.2, run_time=1.6e5, start_time=-1.6e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.01, run_time=6.0e5, start_time=-6.0e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.03, run_time=4.2e5, start_time=-4.2e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.05, run_time=3.5e5, start_time=-3.5e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.1, run_time=2.8e5, start_time=-2.8e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.2, run_time=2.2e5, start_time=-2.2e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.01, run_time=7.5e5, start_time=-7.5e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.03, run_time=5.2e5, start_time=-5.2e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.05, run_time=4.4e5, start_time=-4.4e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.1, run_time=3.5e5, start_time=-3.5e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.2, run_time=2.8e5, start_time=-2.8e5),
]
filenames = map(ips) do ip
    make_filename(ip; θτ=θτ_from_str(ip), βα="init")
end
destfilenamess = map(ips) do ip
    [
        make_filename(ip; θτ=θτ_from_str(ip), βα=0.1), 
        make_filename(ip; θτ=θτ_from_str(ip), βα=0.2), 
        make_filename(ip; θτ=θτ_from_str(ip), βα=0.05)
    ]
end

for (ip, filename, copy_to) in zip(ips, filenames, destfilenamess)
    save_script(filename, T, ip, filename; loc=path, copy_to)
end

println()
println("VARYING COOLING, DEPTH AND STRAIN")
path = "jobs-simulation/cooling-depth"
mkpath(path)
# Run with βα = 0.1
ips = map(ips) do ip
    (; ip..., run_time=5e5, βα=0.1, save_time=3e3)
end

# Varying strain rate for two cooling rates
for (T, ip, destfilenames) in zip(Ts, ips, destfilenamess)
    save_script(destfilenames[1], "04:00:00", (; ip..., βα=0.1, run_time=10e5), destfilenames[1]; loc=path)
    save_script(destfilenames[2], "04:00:00", (; ip..., βα=0.2, run_time=5e5), destfilenames[2]; loc=path)
    save_script(destfilenames[3], "04:00:00", (; ip..., βα=0.05, run_time=20e5), destfilenames[3]; loc=path)
end

# Varying along-front winds

# Varying across-front winds
