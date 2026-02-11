#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=3:00:00
#SBATCH --job-name=Init-b1_0-C1_8
#SBATCH --output=../scratch/logs/Init-b1_0-C1_8.txt

module load julia/1.10.10

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-trig
cd ~/turbulence-at-many-fronts

output_folder=$SCRATCH/turbulence-at-many-fronts/Init-b1_0-C1_8

run_time="4.4e5"
start_time="-4.4e5"
save_time="1e4"

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

Roalpha=0
C=1.8
beta0=8

comment="Initialized front with C = $C, beta b = $betab"

julia -t 24 -- src-simulation/simulation.jl $output_folder $run_time $start_time $save_time $f $L $betax $betah $Nx $Nh $Ny $Nz $betab $betal $betaH $Roalpha $C $beta0 $comment

# Copy to a t > 0 run
mkdir $output_folder/../a0_10-b1_0-C1_8
cp $output_folder/checkpoint_* $output_folder/../a0_10-b1_0-C1_8
