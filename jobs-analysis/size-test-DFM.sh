#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=0:30:00
#SBATCH --job-name=size-test-DFM
#SBATCH --output=../scratch/logs/size-test-DFM.txt

module load julia/1.12.5

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

PPFILE=DFM
RAM=/dev/shm/size-test-DFM
mkdir $RAM
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/size-test-h2/OUTPUT.jld2 $PPFILE DFM.jld2 $RAM/size-test-h2 &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/size-test-h3/OUTPUT.jld2 $PPFILE DFM.jld2 $RAM/size-test-h3 &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/size-test-h4/OUTPUT.jld2 $PPFILE DFM.jld2 $RAM/size-test-h4 &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/size-test-h5/OUTPUT.jld2 $PPFILE DFM.jld2 $RAM/size-test-h5 &
wait

rmdir $RAM
