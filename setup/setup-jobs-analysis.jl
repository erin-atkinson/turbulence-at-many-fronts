using JLD2

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
    outpath = "../scratch/turbulence-at-many-fronts/$foldername/$outputfilename"
    if isfile(outpath)
        iterations, times, sp = jldopen(outpath) do file
            iterations = keys(file["timeseries/t"])
            times = [file["timeseries/t/$iteration"] for iteration in iterations]
            sp = file["metadata/parameters"]
            iterations, times, sp
        end
        println("$foldername: $filename -> $outputfilename ($scriptname) done up to ft = $(times[end] * sp.f)")
    end
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
    println(jobname)
    write(joinpath(loc, jobname * ".sh"), make_script(jobname, foldernames, scriptname, T; filename, outputfilename, job))
    println()
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
    "cooling-only",
    "wind",
    "high-resolution"
]

sets = [
    cooling_depth_01,
    cooling_depth_02,
    cooling_depth_005,
    cooling_only,
    wind,
    highresolution,
]

let set = test_set,
    setname = "test-set"
    save_script("$setname-MEAN", set.filenames, "MEAN", "0:00:00"; loc="jobs-analysis", filename="TIMEAVERAGE.jld2", job=false)
    save_script("$setname-GRADIENTS", set.filenames, "GRADIENTS", "0:00:00"; loc="jobs-analysis", filename="MEAN.jld2", job=false)
    save_script("$setname-ENERGY", set.filenames, "ENERGY", "0:00:00"; loc="jobs-analysis", filename="MEAN.jld2", job=false)
    save_script("$setname-PV", set.filenames, "PV", "0:00:00"; loc="jobs-analysis", filename="MEAN.jld2", job=false)
    save_script("$setname-UBALANCE", set.filenames, "UBALANCE", "2:30:00"; loc="jobs-analysis", filename="MEAN.jld2", job=false)
    save_script("$setname-VBALANCE", set.filenames, "VBALANCE", "2:30:00"; loc="jobs-analysis", filename="MEAN.jld2", job=false)
    save_script("$setname-BBALANCE", set.filenames, "BBALANCE", "2:30:00"; loc="jobs-analysis", filename="MEAN.jld2", job=false)
    save_script("$setname-SLICES", set.filenames, "SLICES", "2:30:00"; loc="jobs-analysis", filename="OUTPUT.jld2", job=false)
    
    save_script("$setname-MEAN-20", set.filenames, "MEAN", "0:00:00"; loc="jobs-analysis", filename="MEAN.jld2", job=false, outputfilename="MEAN-20.jld2")
    save_script("$setname-GRADIENTS-20", set.filenames, "GRADIENTS", "0:00:00"; loc="jobs-analysis", filename="MEAN-20.jld2", job=false, outputfilename="GRADIENTS-20.jld2")
    save_script("$setname-ENERGY-20", set.filenames, "ENERGY", "0:00:00"; loc="jobs-analysis", filename="MEAN-20.jld2", job=false, outputfilename="ENERGY-20.jld2")
    save_script("$setname-PV-20", set.filenames, "PV", "0:00:00"; loc="jobs-analysis", filename="MEAN-20.jld2", job=false, outputfilename="PV-20.jld2")
    save_script("$setname-UBALANCE-20", set.filenames, "UBALANCE", "2:30:00"; loc="jobs-analysis", filename="MEAN-20.jld2", job=false, outputfilename="UBALANCE-20.jld2")
    save_script("$setname-VBALANCE-20", set.filenames, "VBALANCE", "2:30:00"; loc="jobs-analysis", filename="MEAN-20.jld2", job=false, outputfilename="VBALANCE-20.jld2")
    save_script("$setname-BBALANCE-20", set.filenames, "BBALANCE", "2:30:00"; loc="jobs-analysis", filename="MEAN-20.jld2", job=false, outputfilename="BBALANCE-20.jld2")
end

for (set, setname) in zip(sets, setnames)
    save_script("$setname-MEAN", set.filenames, "MEAN", "1:00:00"; loc="jobs-analysis", filename="TIMEAVERAGE.jld2", job=false)
    save_script("$setname-GRADIENTS", set.filenames, "GRADIENTS", "1:00:00"; loc="jobs-analysis", filename="MEAN.jld2", job=false)
    save_script("$setname-ENERGY", set.filenames, "ENERGY", "2:00:00"; loc="jobs-analysis", filename="MEAN.jld2", job=false)
    save_script("$setname-PV", set.filenames, "PV", "2:00:00"; loc="jobs-analysis", filename="MEAN.jld2", job=false)
    save_script("$setname-UBALANCE", set.filenames, "UBALANCE", "2:30:00"; loc="jobs-analysis", filename="MEAN.jld2", job=false)
    save_script("$setname-VBALANCE", set.filenames, "VBALANCE", "2:30:00"; loc="jobs-analysis", filename="MEAN.jld2", job=false)
    save_script("$setname-BBALANCE", set.filenames, "BBALANCE", "2:30:00"; loc="jobs-analysis", filename="MEAN.jld2", job=false)
    save_script("$setname-SLICES", set.filenames, "SLICES", "2:30:00"; loc="jobs-analysis", filename="OUTPUT.jld2", job=false)
    
    save_script("$setname-MEAN-20", set.filenames, "MEAN", "1:00:00"; loc="jobs-analysis", filename="MEAN.jld2", outputfilename="MEAN-20.jld2", job=false)
    save_script("$setname-GRADIENTS-20", set.filenames, "GRADIENTS", "1:00:00"; loc="jobs-analysis", filename="MEAN-20.jld2", outputfilename="GRADIENTS-20.jld2", job=false)
    save_script("$setname-ENERGY-20", set.filenames, "ENERGY", "1:00:00"; loc="jobs-analysis", filename="MEAN-20.jld2", outputfilename="ENERGY-20.jld2", job=false)
    save_script("$setname-PV-20", set.filenames, "PV", "1:00:00"; loc="jobs-analysis", filename="MEAN-20.jld2", outputfilename="PV-20.jld2", job=false)
    save_script("$setname-UBALANCE-20", set.filenames, "UBALANCE", "2:30:00"; loc="jobs-analysis", filename="MEAN-20.jld2", outputfilename="UBALANCE-20.jld2", job=false)
    save_script("$setname-VBALANCE-20", set.filenames, "VBALANCE", "2:30:00"; loc="jobs-analysis", filename="MEAN-20.jld2", outputfilename="VBALANCE-20.jld2", job=false)
    save_script("$setname-BBALANCE-20", set.filenames, "BBALANCE", "2:30:00"; loc="jobs-analysis", filename="MEAN-20.jld2", outputfilename="BBALANCE-20.jld2", job=false)
end

let set = highresolution_init,
    setname = "high-resolution-init"
    save_script("$setname-MEAN", set.filenames, "MEAN", "0:00:00"; loc="jobs-analysis", filename="TIMEAVERAGE.jld2", job=false)
    save_script("$setname-GRADIENTS", set.filenames, "GRADIENTS", "0:00:00"; loc="jobs-analysis", filename="TIMEAVERAGE.jld2", job=false)
    save_script("$setname-ENERGY", set.filenames, "ENERGY", "0:00:00"; loc="jobs-analysis", filename="TIMEAVERAGE.jld2", job=false)
    save_script("$setname-PV", set.filenames, "PV", "0:00:00"; loc="jobs-analysis", filename="TIMEAVERAGE.jld2", job=false)
end
