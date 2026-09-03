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
    """
end

function make_body(foldername, N_window; job=true)
    endstr = job ? "&" : ""
    """
    julia -t 24 -- src-analysis/postprocess/increase_window.jl \$SCRATCH/turbulence-at-many-fronts/$foldername/AVERAGE.jld2 $N_window $endstr
    """
end

function make_cleanup()
    """
    wait
    """
end

function make_script(jobname, foldernames, N_window, T; job=true)
    body = mapreduce(*, foldernames) do foldername
        make_body(foldername, N_window; job)
    end
    return make_preamble(jobname, T) * body * make_cleanup()
end

function save_script(jobname, foldernames, N_window, T; loc="", job=true)
    write(joinpath(loc, jobname * ".sh"), make_script(jobname, foldernames, N_window, T; job))
    return nothing
end

include("ensemble.jl")

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
]

let set = test_set,
    setname = "test-set"
    save_script("$setname-WINDOW", set.filenames, 5, "0:30:00"; loc="jobs-analysis", job=false)
    save_script("$setname-WINDOWWIDE", set.filenames, 20, "0:30:00"; loc="jobs-analysis", job=false)
end

#=
map(sets, setnames) do set, setname
    save_script("$setname-WINDOW", set.filenames, 5, "0:30:00"; loc="jobs-analysis")
    save_script("$setname-WINDOWWIDE", set.filenames, 50, "0:30:00"; loc="jobs-analysis")
end
=#