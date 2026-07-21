#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=0:30:00
#SBATCH --job-name=test-set-increase-window
#SBATCH --output=../scratch/logs/test-set-increase-window.txt

module load julia/1.12.5

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

RAM=/dev/shm/test-set-increase-window
mkdir $RAM
julia -t 24 -- src-analysis/postprocess/increase_window.jl $SCRATCH/turbulence-at-many-fronts//OUTPUT.jld2 3 $RAM/size-test-h2 &
wait

rm $RAM -rf
