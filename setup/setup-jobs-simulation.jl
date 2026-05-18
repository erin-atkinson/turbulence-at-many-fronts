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

    julia -t 8 -- src-simulation/simulation.jl \$output_folder \$run_time \$start_time \$max_time \$save_time \$f \$L \$betax \$betah \$Nx \$Nh \$Ny \$Nz \$betal \$betaH \$betaalpha \$betaB \$betatau \$thetatau \$beta0 \$comment &
    
    """
end

function make_cleanup()
    """
    wait
    """
end

function make_copy(filename::String)
    """
    
    mkdir \$output_folder/../$filename
    cp \$output_folder/checkpoint* \$output_folder/../$filename
    """
end

function make_copy(filenames)
    mapreduce(make_copy, *, filenames)
end


make_script(jobname, T, ip::NamedTuple, filename::String) = make_script(jobname, T, [ip], [filename])

function make_script(jobname, T, ips, filenames, destfilenames=[])
    body = mapreduce(make_body, *, ips, filenames)
    return make_preamble(jobname, T) * body * make_cleanup() * make_copy(destfilenames)
end

function save_script(jobname, T, ips, filenames; loc="", copy_to=[])
    scriptpath = joinpath(loc, jobname * ".sh")
    write(scriptpath, make_script(jobname, T, ips, filenames, move_to))
    return nothing
end

default_ip = (;
    run_time = 8e5, start_time = -2e5, save_time = 3e3, max_time = -1.0,
    f = 1e-4, L = 1e3,
    βx = 20, βh = 4,
    Nx = 1024, Nh = 768, Ny = 128, Nz = 64,
    βℓ = 6, βH = 0.1,
    βα = 0.1, βB = 0.03, βτ = 0.00, θτ = 0.0, β₀ = 6,
    comment = ""
)

# Tests: varying central region size
jobname = "size-test"
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
save_script(jobname, "5:00:00", ips, filenames; loc="jobs-simulation")

# Tests: varying resolution
jobname = "resolution-test"
ips = [
    (; default_ip..., save_time=1e4, Nx=1024, Nh=768, Nz=96, comment="Resolution test x1024 z96"),
    (; default_ip..., save_time=1e4, Nx=1024, Nh=768, Nz=96, comment="Resolution test x1024 z32"),
    (; default_ip..., save_time=1e4, Nx=800, Nh=600, Nz=64, comment="Resolution test x800 z64"),
    (; default_ip..., save_time=1e4, Nx=800, Nh=600, Nz=32, comment="Resolution test x800 z32"),
]
filenames = [
    "resolution-test-x1024-z96",
    "resolution-test-x1024-z32",
    "resolution-test-x800-z64",
    "resolution-test-x800-z32",
]
save_script(jobname, "5:00:00", ips, filenames; loc="jobs-simulation")

# Varying cooling rate and depth
# Initialise
# Run with βα = 0.1

# Varying strain rate for two cooling rates

# Varying along-front winds

# Varying across-front winds
