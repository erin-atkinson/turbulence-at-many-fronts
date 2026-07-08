#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=4:00:00
#SBATCH --job-name=outer-test-x900
#SBATCH --output=../scratch/logs/outer-test-x900.txt

module load julia/1.12.5

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-trig
# export JULIA_CUDA_SOFT_MEMORY_LIMIT=10%
cd ~/turbulence-at-many-fronts

output_folder=$SCRATCH/turbulence-at-many-fronts/outer-test-x900

run_time=800000.0
start_time=-200000.0
max_time=-1.0
save_time=10000.0

f=0.0001
L=1000.0

betax=20
betah=3

Nx=900
Nh=700
Ny=128
Nz=64

betal=6
betaH=0.1

betaalpha=0.1
betaB=0.03
betatau=0.0
thetatau=0.0
beta0=6

comment="Outer test x900"

julia -t 8 -- src-simulation/simulation.jl $output_folder $run_time $start_time $max_time $save_time $f $L $betax $betah $Nx $Nh $Ny $Nz $betal $betaH $betaalpha $betaB $betatau $thetatau $beta0 $comment
