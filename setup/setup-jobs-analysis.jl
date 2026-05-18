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

foldernames = [
    "size-test-h2",
    "size-test-h3",
    "size-test-h4",
    "size-test-h5"
]

save_script("size-test-DFM", foldernames, "DFM", "0:30:00"; loc="jobs-analysis")
