#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=0:30:00
#SBATCH --job-name=ppDFM
#SBATCH --output=../scratch/logs/ppDFM.txt

module load julia/1.10.10

export RAM=/dev/shm/DFM
mkdir $RAM

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/Strain
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER DFM $RAM

rm $RAM -rf
