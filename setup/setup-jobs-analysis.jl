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

function make_body(foldername, scriptname; filename="OUTPUT.jld2", outputfilename="$scriptname.jld2")
    """
    julia -t 24 -- src-analysis/postprocess/postprocess.jl \$SCRATCH/turbulence-at-many-fronts/$foldername/$filename \$PPFILE $outputfilename \$RAM/$foldername &
    """
end

function make_cleanup()
    """
    wait
    
    rmdir \$RAM
    """
end

function make_script(jobname, foldernames, scriptname, T; filename="OUTPUT.jld2", outputfilename="$scriptname.jld2")
    body = mapreduce(*, foldernames) do foldername
        make_body(foldername, scriptname; filename, outputfilename)
    end
    return make_preamble(jobname, scriptname, T) * body * make_cleanup()
end

function save_script(jobname, foldernames, scriptname, T; loc="", filename="OUTPUT.jld2", outputfilename="$scriptname.jld2")
    write(joinpath(loc, jobname * ".sh"), make_script(jobname, foldernames, scriptname, T; filename, outputfilename))
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
    "test_set",
    "test_set_init",
]

sets = [
    cooling_depth_01,
    cooling_depth_02,
    cooling_depth_005,
    cooling_only,
    depth_only,
    test_set,
    test_set_init
]

map(sets, setnames) do set, setname
    save_script("$setname-MEAN", set.filenames, "MEAN", "0:30:00"; loc="jobs-analysis")
    save_script("$setname-MEAN-5", set.filenames, "MEAN", "0:30:00"; loc="jobs-analysis", filename="OUTPUT-5.jld2", outputfilename="MEAN-5.jld2")
    save_script("$setname-GRADIENTS", set.filenames, "GRADIENTS", "0:30:00"; loc="jobs-analysis")
    
    save_script("$setname-UBALANCE", set.filenames, "UBALANCE", "0:45:00"; loc="jobs-analysis")
    save_script("$setname-VBALANCE", set.filenames, "VBALANCE", "0:45:00"; loc="jobs-analysis")
    save_script("$setname-BBALANCE", set.filenames, "BBALANCE", "0:45:00"; loc="jobs-analysis")
    
    save_script("$setname-UBALANCE-5", set.filenames, "UBALANCE", "0:45:00"; loc="jobs-analysis", filename="OUTPUT-5.jld2", outputfilename="UBALANCE-5.jld2")
    save_script("$setname-VBALANCE-5", set.filenames, "VBALANCE", "0:45:00"; loc="jobs-analysis", filename="OUTPUT-5.jld2", outputfilename="VBALANCE-5.jld2")
    save_script("$setname-BBALANCE-5", set.filenames, "BBALANCE", "0:45:00"; loc="jobs-analysis", filename="OUTPUT-5.jld2", outputfilename="BBALANCE-5.jld2")

    save_script("$setname-ENERGY", set.filenames, "ENERGY", "2:30:00"; loc="jobs-analysis")
    save_script("$setname-ENERGY-5", set.filenames, "ENERGY", "2:30:00"; loc="jobs-analysis", filename="OUTPUT-5.jld2", outputfilename="ENERGY-5.jld2")

    #save_script("$setname-PV", set.filenames, "PV", "1:00:00"; loc="jobs-analysis")
    #save_script("$setname-SLICES", set.filenames, "SLICES", "0:45:00"; loc="jobs-analysis")
    #save_script("$setname-BBALANCE", set.filenames, "BBALANCE", "0:45:00"; loc="jobs-analysis")
end
