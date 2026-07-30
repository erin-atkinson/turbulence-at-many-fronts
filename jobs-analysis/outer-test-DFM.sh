#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=0:30:00
#SBATCH --job-name=outer-test-DFM
#SBATCH --output=../scratch/logs/outer-test-DFM.txt

module load julia/1.12.5

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

PPFILE=DFM
RAM=/dev/shm/outer-test-DFM
mkdir $RAM
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/outer-test-x800/OUTPUT.jld2 $PPFILE DFM.jld2 $RAM/outer-test-x800 &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/outer-test-x900/OUTPUT.jld2 $PPFILE DFM.jld2 $RAM/outer-test-x900 &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/outer-test-x1000/OUTPUT.jld2 $PPFILE DFM.jld2 $RAM/outer-test-x1000 &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/outer-test-x1100/OUTPUT.jld2 $PPFILE DFM.jld2 $RAM/outer-test-x1100 &
wait

rmdir $RAM
