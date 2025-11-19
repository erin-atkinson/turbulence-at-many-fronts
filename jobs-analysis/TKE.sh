#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=4:00:00
#SBATCH --job-name=ppTKE
#SBATCH --output=../scratch/logs/ppTKE.txt

module load julia/1.10.10

# Copy installation to RAM disk
echo "Copying installation to RAM disk"
export RAM=/dev/shm/TKE
mkdir $RAM

cp -r $HOME/turbulence-at-many-fronts/.julia-tri $RAM

# Launch from RAM disk
echo "Running..."
export JULIA_DEPOT_PATH=$RAM/.julia-tri
export JULIA_SCRATCH_TRACK_ACCESS=0
cd ~/turbulence-at-many-fronts

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/Strain
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER TKE $RAM

rm $RAM -rf
