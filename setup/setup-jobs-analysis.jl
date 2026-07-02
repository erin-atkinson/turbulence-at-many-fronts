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

function make_preamble(jobname, scriptname, T)
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
    
    PPFILE=$scriptname
    RAM=/dev/shm/$jobname
    mkdir \$RAM
    """
end

function make_body(foldername, scriptname)
    """
    julia -t 24 -- src-analysis/postprocess/postprocess.jl \$SCRATCH/turbulence-at-many-fronts/$foldername \$PPFILE \$RAM/$foldername &
    """
end

function make_cleanup()
    """
    wait
    
    rm \$RAM -rf
    """
end

function make_script(jobname, foldernames, scriptname, T)
    body = mapreduce(*, foldernames) do foldername
        make_body(foldername, scriptname)
    end
    return make_preamble(jobname, scriptname, T) * body * make_cleanup()
end

function save_script(jobname, foldernames, scriptname, T; loc="")
    write(joinpath(loc, jobname * ".sh"), make_script(jobname, foldernames, scriptname, T))
    return nothing
end

include("ensemble.jl")

save_script("size-test-DFM", size_test.filenames, "DFM", "0:30:00"; loc="jobs-analysis")
save_script("resolution-test-DFM", resolution_test.filenames, "DFM", "0:30:00"; loc="jobs-analysis")
save_script("outer-test-DFM", outer_test.filenames, "DFM", "0:30:00"; loc="jobs-analysis")

save_script("cooling-depth-init-DFM", cooling_depth_init.filenames, "DFM", "0:45:00"; loc="jobs-analysis")
save_script("cooling-depth-init-TKE", cooling_depth_init.filenames, "TKE", "0:45:00"; loc="jobs-analysis")
save_script("cooling-depth-init-TKE-L", cooling_depth_init.filenames, "TKE-L", "0:45:00"; loc="jobs-analysis")

setnames = [
    "cooling-depth-0_1",
    "cooling-depth-0_2",
    "cooling-depth-0_05",
    "cooling-only",
    "depth-only",
]

sets = [
    cooling_depth_01,
    cooling_depth_02,
    cooling_depth_005,
    cooling_only,
    depth_only,
]

map(sets, setnames) do set, setname
    save_script("$setname-DFM", set.filenames, "DFM", "0:45:00"; loc="jobs-analysis")
    save_script("$setname-COARSE", set.filenames, "COARSE", "0:45:00"; loc="jobs-analysis")
    save_script("$setname-TKE", set.filenames, "TKE", "2:30:00"; loc="jobs-analysis")
    save_script("$setname-TKE-L", set.filenames, "TKE-L", "2:30:00"; loc="jobs-analysis")
    save_script("$setname-PV", set.filenames, "PV", "1:00:00"; loc="jobs-analysis")
    save_script("$setname-PV-L", set.filenames, "PV-L", "2:30:00"; loc="jobs-analysis")
    save_script("$setname-SLICES", set.filenames, "SLICES", "0:45:00"; loc="jobs-analysis")
end
