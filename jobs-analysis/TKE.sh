#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=4:00:00
#SBATCH --job-name=ppTKE
#SBATCH --output=../scratch/logs/ppTKE.txt

module load julia/1.10.10

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/Strain
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER TKE $RAM

rm $RAM -rf
