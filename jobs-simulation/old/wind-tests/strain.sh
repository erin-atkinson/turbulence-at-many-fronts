#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=3:00:00
#SBATCH --job-name=wind-test-strain
#SBATCH --output=../scratch/logs/wind-test-strain.txt

module load julia/1.10.10

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-trig
cd ~/turbulence-at-many-fronts

output_folder=$SCRATCH/turbulence-at-many-fronts/wind-test-strain

run_time="8e5"
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
C=0.0
W=-0.07
beta0=8

comment="Some windy strain"

julia -t 24 -- src-simulation/simulation.jl $output_folder $run_time $start_time $save_time $f $L $betax $betah $Nx $Nh $Ny $Nz $betab $betal $betaH $Roalpha $C $W $beta0 $comment
