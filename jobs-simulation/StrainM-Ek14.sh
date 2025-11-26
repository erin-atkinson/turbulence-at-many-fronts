#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=24:00:00
#SBATCH --job-name=StrainM-Ek14
#SBATCH --output=../scratch/logs/Strain-%j.txt

module load julia/1.10.10

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-trig
cd ~/turbulence-at-many-fronts

output_folder=$SCRATCH/turbulence-at-many-fronts/StrainM-Ek14

run_time="10e5"
start_time="-6e5"
save_time="1e3"

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

Roalpha=0.1
Ek=1.4
beta0=8

comment="Frontogenesis with alpha = 0.1f, Ek = 1.4"

julia -t 24 -- src-simulation/simulation.jl $output_folder $run_time $start_time $save_time $f $L $betax $betah $Nx $Nh $Ny $Nz $betab $betal $betaH $Roalpha $Ek $beta0 $comment
