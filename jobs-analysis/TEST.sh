#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=0:30:00
#SBATCH --job-name=cooling-only-TEST
#SBATCH --output=../scratch/logs/cooling-only-TEST.txt

module load julia/1.12.5

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

PPFILE=TEST
RAM=/dev/shm/cooling-only-TEST
mkdir $RAM

julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_05-0_0-C $PPFILE $RAM/0_1-0_1-0_05-0_0-C &
wait

rm $RAM -rf
