#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=2:15:00
#SBATCH --job-name=ppPV
#SBATCH --output=../scratch/logs/ppPV.txt

module load julia/1.10.10

# Copy installation to RAM disk
echo "Copying installation to RAM disk"
export RAM=/dev/shm/PV
mkdir $RAM

cp -r $HOME/turbulence-at-many-fronts/.julia-tri $RAM

# Launch from RAM disk
echo "Running..."
export JULIA_DEPOT_PATH=$RAM/.julia-tri
export JULIA_SCRATCH_TRACK_ACCESS=0
cd ~/turbulence-at-many-fronts

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/Strain
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER PV $RAM

rm $RAM -rf
