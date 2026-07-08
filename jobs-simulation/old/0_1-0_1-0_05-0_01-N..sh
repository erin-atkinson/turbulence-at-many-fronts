#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=00:00:00
#SBATCH --job-name=0_1-0_1-0_05-0_01-N
#SBATCH --output=../scratch/logs/0_1-0_1-0_05-0_01-N.txt

module load julia/1.12.5

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-trig
cd ~/turbulence-at-many-fronts

output_folder=$SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_05-0_01-N

run_time=800000.0
start_time=-400000.0
save_time=1000.0

f=0.0001
L=1000.0

betax=20
betah=3

Nx=1024
Nh=768
Ny=128
Nz=64

betal=6
betaH=0.1

betaalpha=0.1
betaB=0.05
betatau=0.01
betatheta=0.0
beta0=8

comment=

julia -t 24 -- src-simulation/simulation.jl $output_folder $run_time $start_time $save_time $f $L $betax $betah $Nx $Nh $Ny $Nz $betal $betaH $betaalpha $betaB $betatau $thetatau $beta0 $comment
