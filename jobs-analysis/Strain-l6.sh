#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=4:00:00
#SBATCH --job-name=ppStrain-l6
#SBATCH --output=../scratch/logs/ppStrain-l6.txt

module load julia/1.10.10

# Copy installation to RAM disk
echo "Copying installation to RAM disk"
export RAM=/dev/shm/Strain-l6
mkdir $RAM

cp -r $HOME/turbulence-at-many-fronts/.julia-tri $RAM

# Launch from RAM disk
echo "Running..."
export JULIA_DEPOT_PATH=$RAM/.julia-tri
export JULIA_SCRATCH_TRACK_ACCESS=0
cd ~/turbulence-at-many-fronts

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-many-fronts/Strain-l6
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER DFM $RAM
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER PV $RAM
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER RI $RAM
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER TKE $RAM

rm $RAM -rf
