#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=6:00:00
#SBATCH --job-name=Strain-l8
#SBATCH --output=../scratch/logs/Strain-%j.txt

module load julia/1.10.10

# Copy installation to RAM disk
export RAM=/dev/shm/Strain-l8
mkdir $RAM

cp -rn $HOME/turbulence-at-many-fronts/.julia-trig $RAM

# Launch from RAM disk
export JULIA_DEPOT_PATH=$RAM/.julia-trig
export JULIA_SCRATCH_TRACK_ACCESS=0
cd ~/turbulence-at-many-fronts

output_folder=$SCRATCH/turbulence-at-many-fronts/Strain-l8

run_time="16e5"
start_time="-4e5"
save_time="2e3"

f="1e-4"
L="1e3"

betax=32
betah=4

Nx=1024
Nh=768
Ny=128
Nz=64

betab=1
betal=8
betaH=0.1

alpha="1e-5"
Q=100
beta0=8

comment="Checking Ro dependence"

julia -t 24 -- src-simulation/simulation.jl $output_folder $run_time $start_time $save_time $f $L $betax $betah $Nx $Nh $Ny $Nz $betab $betal $betaH $alpha $Q $beta0 $comment
