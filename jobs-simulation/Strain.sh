#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=24:00:00
#SBATCH --job-name=Strain
#SBATCH --output=../scratch/logs/Strain-%j.txt

module load julia/1.10.10

# Copy installation to RAM disk
export RAM=/dev/shm/turbulence-at-many-fronts
mkdir $RAM

cp -r $HOME/turbulence-at-many-fronts/.julia-trig $RAM

# Launch from RAM disk
export JULIA_DEPOT_PATH=$RAM/.julia-trig
export JULIA_SCRATCH_TRACK_ACCESS=0
cd ~/turbulence-at-many-fronts

output_folder=$SCRATCH/turbulence-at-many-fronts/Strain

run_time="4e5"
start_time="-4e5"
save_time="1e3"

f="1e-4"

Lx=4
Lh=0.5
H=100

Nx=1024
Nh=768
Ny=128
Nz=128

Ro=0.1
Ri=2.0
db=0.029

alpha="1e-5"
Q=100
N0sq=0.003

comment="Checking Ro dependence"

julia -t 24 -- src-simulation/simulation.jl $output_folder $run_time $start_time $save_time $f $Lx $Lh $H $Nx $Nh $Ny $Nz $Ro $Ri $db $alpha $Q $N0sq $comment
