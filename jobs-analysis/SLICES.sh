#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=1:00:00
#SBATCH --job-name=ppSLICES
#SBATCH --output=../scratch/logs/ppSLICES.txt

module load julia/1.10.10

export RAM=/dev/shm/SLICES
mkdir $RAM

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/Strain-l8
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER SLICES $RAM

rm $RAM -rf
