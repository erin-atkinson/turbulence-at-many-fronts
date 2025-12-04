#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=5:00:00
#SBATCH --job-name=ppSCALES
#SBATCH --output=../scratch/logs/ppSCALES.txt

module load julia/1.10.10

export RAM=/dev/shm/SCALES
mkdir $RAM

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/StrainM-C07
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER SCALES $RAM

export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/StrainM-C10
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER SCALES $RAM

export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/StrainM-C14
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER SCALES $RAM

export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/StrainM-C20
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER SCALES $RAM

rm $RAM -rf
