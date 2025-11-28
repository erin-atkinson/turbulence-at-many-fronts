#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=4:00:00
#SBATCH --job-name=Init-C07
#SBATCH --output=../scratch/logs/Strain-%j.txt

module load julia/1.10.10

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-trig
cd ~/turbulence-at-many-fronts

output_folder=$SCRATCH/turbulence-at-many-fronts/Init-C07

run_time="12e5"
start_time="-12e5"
save_time="1e4"

f="1e-4"
L="1e3"

betax=20
betah=3

Nx=1024
Nh=768
Ny=128
Nz=128

betab=1
betal=6
betaH=0.1

Roalpha=0
C=0.7
beta0=8

comment="Initialized front with C = 0.7"

julia -t 24 -- src-simulation/simulation.jl $output_folder $run_time $start_time $save_time $f $L $betax $betah $Nx $Nh $Ny $Nz $betab $betal $betaH $Roalpha $C $beta0 $comment

# Copy to each t > 0 run
mkdir $output_folder/../StrainL-C07
cp $output_folder/checkpoint_* $output_folder/../StrainL-C07

mkdir $output_folder/../StrainM-C07
cp $output_folder/checkpoint_* $output_folder/../StrainM-C07

mkdir $output_folder/../StrainH-C07
cp $output_folder/checkpoint_* $output_folder/../StrainH-C07
