#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=3:00:00
#SBATCH --job-name=ppBBALANCE
#SBATCH --output=../scratch/logs/ppBBALANCE.txt

module load julia/1.10.10

export RAM=/dev/shm/DFM
mkdir $RAM

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/a0_10-b1_0-C1_0
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER BBALANCE $RAM

export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/a0_10-b1_0-C2_0
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER BBALANCE $RAM

export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/a0_10-b1_0-C3_0
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER BBALANCE $RAM

export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/a0_10-b1_0-C1_4
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER BBALANCE $RAM

export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/a0_10-b1_0-C1_8
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER BBALANCE $RAM

#export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/a0_10-b1_0-C2_2
#julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER BBALANCE $RAM

export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/a0_10-b1_0-C2_6
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER BBALANCE $RAM

rm $RAM -rf
