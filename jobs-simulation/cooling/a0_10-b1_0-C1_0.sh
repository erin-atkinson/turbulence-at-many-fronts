#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=6:00:00
#SBATCH --job-name=a0_10-b1_0-C1_0
#SBATCH --output=../scratch/logs/a0_10-b1_0-C1_0.txt

module load julia/1.10.10

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-trig
cd ~/turbulence-at-many-fronts

output_folder=$SCRATCH/turbulence-at-many-fronts/a0_10-b1_0-C1_0

run_time="10e5"
start_time="-8e5"
save_time="1e3"

f="1e-4"
L="1e3"

betax=20
betah=3

Nx=1024
Nh=768
Ny=128
Nz=64

betab=1
betal=6
betaH=0.1

Roalpha=0.1
C=1.0
beta0=8

comment="Frontogenesis with Roalpha = $Roalpha, C = $C, beta b = $betab"

julia -t 24 -- src-simulation/simulation.jl $output_folder $run_time $start_time $save_time $f $L $betax $betah $Nx $Nh $Ny $Nz $betab $betal $betaH $Roalpha $C $beta0 $comment
