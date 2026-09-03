
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

function make_body(foldername, scriptname; filename="OUTPUT.jld2", outputfilename="$scriptname.jld2", job=true)
    threadstr = job ? "24" : "2"
    ramstr = job ? "\$RAM/$foldername &" : ""
    """
    julia -t $threadstr -- src-analysis/postprocess/postprocess.jl \$SCRATCH/turbulence-at-many-fronts/$foldername/$filename \$PPFILE $outputfilename $ramstr
    """
end

function make_cleanup()
    """
    wait
    
    rmdir \$RAM
    """
end

function make_script(jobname, foldernames, scriptname, T; filename="OUTPUT.jld2", outputfilename="$scriptname.jld2", job=true)
    body = mapreduce(*, foldernames) do foldername
        make_body(foldername, scriptname; filename, outputfilename, job)
    end
    return make_preamble(jobname, scriptname, T) * body * make_cleanup()
end

function save_script(jobname, foldernames, scriptname, T; loc="", filename="OUTPUT.jld2", outputfilename="$scriptname.jld2", job=true)
    write(joinpath(loc, jobname * ".sh"), make_script(jobname, foldernames, scriptname, T; filename, outputfilename, job))
    return nothing
end

include("ensemble.jl")

save_script("size-test-DFM", size_test.filenames, "DFM", "0:30:00"; loc="jobs-analysis")
save_script("resolution-test-DFM", resolution_test.filenames, "DFM", "0:30:00"; loc="jobs-analysis")
save_script("outer-test-DFM", outer_test.filenames, "DFM", "0:30:00"; loc="jobs-analysis")

setnames = [
    "cooling-depth-0_1",
    "cooling-depth-0_2",
    "cooling-depth-0_05",
    "wind"
]

sets = [
    cooling_depth_01,
    cooling_depth_02,
    cooling_depth_005,
    wind
]

let set = test_set,
    setname = "test-set"
    save_script("$setname-MEAN", set.filenames, "MEAN", "0:00:00"; loc="jobs-analysis", filename="AVERAGE.jld2", job=false)
    save_script("$setname-GRADIENTS", set.filenames, "GRADIENTS", "0:00:00"; loc="jobs-analysis", filename="AVERAGE.jld2", job=false)
    save_script("$setname-ENERGY", set.filenames, "ENERGY", "0:00:00"; loc="jobs-analysis", filename="AVERAGE.jld2", job=false)
    save_script("$setname-PV", set.filenames, "PV", "0:00:00"; loc="jobs-analysis", filename="AVERAGE.jld2", job=false)
    save_script("$setname-STREAMFUNCTION", set.filenames, "STREAMFUNCTION", "2:30:00"; loc="jobs-analysis", filename="AVERAGE.jld2", job=false, outputfilename="STREAMFUNCTION.jld2")
    
    save_script("$setname-MEAN-20", set.filenames, "MEAN", "0:00:00"; loc="jobs-analysis", filename="AVERAGE-20.jld2", job=false, outputfilename="MEAN-20.jld2")
    save_script("$setname-GRADIENTS-20", set.filenames, "GRADIENTS", "0:00:00"; loc="jobs-analysis", filename="AVERAGE-20.jld2", job=false, outputfilename="GRADIENTS-20.jld2")
    save_script("$setname-ENERGY-20", set.filenames, "ENERGY", "0:00:00"; loc="jobs-analysis", filename="AVERAGE-20.jld2", job=false, outputfilename="ENERGY-20.jld2")
    save_script("$setname-PV-20", set.filenames, "PV", "0:00:00"; loc="jobs-analysis", filename="AVERAGE-20.jld2", job=false, outputfilename="PV-20.jld2")
    save_script("$setname-STREAMFUNCTION-20", set.filenames, "STREAMFUNCTION", "0:00:00"; loc="jobs-analysis", filename="AVERAGE-20.jld2", job=false, outputfilename="STREAMFUNCTION-20.jld2")
end

#=
map(sets, setnames) do set, setname
    save_script("$setname-MEAN", set.filenames, "MEAN", "0:30:00"; loc="jobs-analysis", filename="AVERAGE.jld2")
    save_script("$setname-GRADIENTS", set.filenames, "GRADIENTS", "0:30:00"; loc="jobs-analysis", filename="AVERAGE.jld2")
    
    save_script("$setname-UBALANCE", set.filenames, "UBALANCE", "0:45:00"; loc="jobs-analysis", filename="AVERAGE.jld2")
    save_script("$setname-VBALANCE", set.filenames, "VBALANCE", "0:45:00"; loc="jobs-analysis", filename="AVERAGE.jld2")
    save_script("$setname-BBALANCE", set.filenames, "BBALANCE", "0:45:00"; loc="jobs-analysis", filename="AVERAGE.jld2")
    
    #save_script("$setname-UBALANCE-5", set.filenames, "UBALANCE", "0:45:00"; loc="jobs-analysis", filename="OUTPUT-5.jld2", outputfilename="UBALANCE-5.jld2")
    #save_script("$setname-VBALANCE-5", set.filenames, "VBALANCE", "0:45:00"; loc="jobs-analysis", filename="OUTPUT-5.jld2", outputfilename="VBALANCE-5.jld2")
    #save_script("$setname-BBALANCE-5", set.filenames, "BBALANCE", "0:45:00"; loc="jobs-analysis", filename="OUTPUT-5.jld2", outputfilename="BBALANCE-5.jld2")

    save_script("$setname-ENERGY", set.filenames, "ENERGY", "2:30:00"; loc="jobs-analysis", filename="AVERAGE.jld2")
    #save_script("$setname-ENERGY-5", set.filenames, "ENERGY", "2:30:00"; loc="jobs-analysis", filename="OUTPUT-5.jld2", outputfilename="ENERGY-5.jld2")

    #save_script("$setname-GRADIENTS-5", set.filenames, "GRADIENTS", "0:30:00"; loc="jobs-analysis", filename="OUTPUT-5.jld2", outputfilename="GRADIENTS-5.jld2")

    #save_script("$setname-SURFACE", set.filenames, "SURFACE", "0:30:00"; loc="jobs-analysis", filename="OUTPUT-20.jld2")
    
    #save_script("$setname-PV", set.filenames, "PV", "1:00:00"; loc="jobs-analysis")
    #save_script("$setname-SLICES", set.filenames, "SLICES", "0:45:00"; loc="jobs-analysis")
    #save_script("$setname-BBALANCE", set.filenames, "BBALANCE", "0:45:00"; loc="jobs-analysis")
end
=#