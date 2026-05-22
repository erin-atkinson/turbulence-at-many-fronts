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

foldernames = [
    "resolution-test-x1400-z96",
    "resolution-test-x1400-z64",
    "resolution-test-x1024-z96",
    "resolution-test-x1024-z32",
    "resolution-test-x800-z64",
    "resolution-test-x800-z32",
]
save_script("resolution-test-DFM", foldernames, "DFM", "0:30:00"; loc="jobs-analysis")

foldernames = [
    "outer-test-x800",
    "outer-test-x900",
    "outer-test-x1000",
    "outer-test-x1100",
]
save_script("outer-test-DFM", foldernames, "DFM", "0:30:00"; loc="jobs-analysis")

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
foldernames = map(ips) do ip
    make_filename(ip; θτ=θτ_from_str(ip), βα="init")
end
save_script("cooling-depth-init-DFM", foldernames, "DFM", "0:45:00"; loc="jobs-analysis")
save_script("cooling-depth-init-TKE", foldernames, "TKE", "0:45:00"; loc="jobs-analysis")
save_script("cooling-depth-init-TKE-L", foldernames, "TKE-L", "0:45:00"; loc="jobs-analysis")

