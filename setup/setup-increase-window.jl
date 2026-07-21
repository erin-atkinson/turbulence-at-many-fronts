function make_filename(sp, ext=nothing, pre=""; βH=sp.βH, βα=sp.βα, βB=sp.βB, βτ=sp.βτ, θτ=sp.θτ)
    strs = map([βH, βα, βB, βτ, θτ]) do β
        replace(string(β), "."=>"_")
    end
    ext = isnothing(ext) ? "" : ".$ext"
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
    #SBATCH --ntasks-per-node=192
    #SBATCH --time=$T
    #SBATCH --job-name=$jobname
    #SBATCH --output=../scratch/logs/$jobname.txt

    module load julia/1.12.5

    # Launch from scratch
    export JULIA_DEPOT_PATH=\$SCRATCH/julia-tri
    cd ~/turbulence-at-many-fronts

    RAM=/dev/shm/$jobname
    mkdir \$RAM
    """
end

function make_body(foldername, N_window)
    """
    julia -t 24 -- src-analysis/postprocess/increase_window.jl \$SCRATCH/turbulence-at-many-fronts/$foldername/OUTPUT.jld2 $N_window \$RAM/$foldername &
    """
end

function make_cleanup()
    """
    wait
    
    rm \$RAM -rf
    """
end

function make_script(jobname, foldernames, N_window, T)
    body = mapreduce(*, foldernames) do foldername
        make_body(foldername, N_window)
    end
    return make_preamble(jobname, T) * body * make_cleanup()
end

function save_script(jobname, foldernames, N_window, T; loc="")
    write(joinpath(loc, jobname * ".sh"), make_script(jobname, foldernames, N_window, T))
    return nothing
end

include("ensemble.jl")

save_script("$setname-WINDOW", set.filenames, 5, "0:30:00"; loc="jobs-analysis")

setnames = [
    "cooling-depth-0_1",
    "cooling-depth-0_2",
    "cooling-depth-0_05",
    "cooling-only",
    "depth-only",
    "test_set"
]

sets = [
    cooling_depth_01,
    cooling_depth_02,
    cooling_depth_005,
    cooling_only,
    depth_only,
    test_set
]

map(sets, setnames) do set, setname
    save_script("$setname-WINDOW", set.filenames, 5, "0:30:00"; loc="jobs-analysis")
end
