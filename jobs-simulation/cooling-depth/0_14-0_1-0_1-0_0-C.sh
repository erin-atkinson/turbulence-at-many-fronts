#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=4:00:00
#SBATCH --job-name=0_14-0_1-0_1-0_0-C
#SBATCH --output=../scratch/logs/0_14-0_1-0_1-0_0-C.txt

module load julia/1.12.5

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-trig
# export JULIA_CUDA_SOFT_MEMORY_LIMIT=10%
cd ~/turbulence-at-many-fronts

output_folder=$SCRATCH/turbulence-at-many-fronts/0_14-0_1-0_1-0_0-C

stop_time=1.0e6
start_time=-350000.0
max_time=-1.0
save_time=10000.0

f=0.0001
L=1000.0

betax=20
betah=3

Nx=1400
Nh=1024
Ny=128
Nz=90

betal=6
betaH=0.14

betaalpha=0.1
betaB=0.1
betatau=0.0
thetatau=0.0
beta0=60

comment="Mixed layer frontogenesis"

julia -t 8 -- src-simulation/simulation.jl $output_folder $stop_time $start_time $max_time $save_time $f $L $betax $betah $Nx $Nh $Ny $Nz $betal $betaH $betaalpha $betaB $betatau $thetatau $beta0 $comment
