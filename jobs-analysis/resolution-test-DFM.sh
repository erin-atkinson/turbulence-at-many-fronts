#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=0:30:00
#SBATCH --job-name=resolution-test-DFM
#SBATCH --output=../scratch/logs/resolution-test-DFM.txt

module load julia/1.12.5

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

PPFILE=DFM
RAM=/dev/shm/resolution-test-DFM
mkdir $RAM
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/resolution-test-x1400-z96 $PPFILE $RAM/resolution-test-x1400-z96 &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/resolution-test-x1400-z64 $PPFILE $RAM/resolution-test-x1400-z64 &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/resolution-test-x1024-z96 $PPFILE $RAM/resolution-test-x1024-z96 &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/resolution-test-x1024-z32 $PPFILE $RAM/resolution-test-x1024-z32 &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/resolution-test-x800-z64 $PPFILE $RAM/resolution-test-x800-z64 &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/resolution-test-x800-z32 $PPFILE $RAM/resolution-test-x800-z32 &
wait

rm $RAM -rf
