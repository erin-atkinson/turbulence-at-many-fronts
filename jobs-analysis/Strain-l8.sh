#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=4:00:00
#SBATCH --job-name=ppStrain-l8
#SBATCH --output=../scratch/logs/ppStrain-l8.txt

module load julia/1.10.10

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/Strain-l8
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER DFM $RAM
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER PV $RAM
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER RI $RAM
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER TKE $RAM

rm $RAM -rf
