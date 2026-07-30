#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=0:45:00
#SBATCH --job-name=depth-only-BBALANCE
#SBATCH --output=../scratch/logs/depth-only-BBALANCE.txt

module load julia/1.12.5

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

PPFILE=BBALANCE
RAM=/dev/shm/depth-only-BBALANCE
mkdir $RAM
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_06-0_1-0_03-0_0-C/OUTPUT.jld2 $PPFILE BBALANCE.jld2 $RAM/0_06-0_1-0_03-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_08-0_1-0_03-0_0-C/OUTPUT.jld2 $PPFILE BBALANCE.jld2 $RAM/0_08-0_1-0_03-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_09-0_1-0_03-0_0-C/OUTPUT.jld2 $PPFILE BBALANCE.jld2 $RAM/0_09-0_1-0_03-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_03-0_0-C/OUTPUT.jld2 $PPFILE BBALANCE.jld2 $RAM/0_1-0_1-0_03-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_12-0_1-0_03-0_0-C/OUTPUT.jld2 $PPFILE BBALANCE.jld2 $RAM/0_12-0_1-0_03-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_13-0_1-0_03-0_0-C/OUTPUT.jld2 $PPFILE BBALANCE.jld2 $RAM/0_13-0_1-0_03-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_14-0_1-0_03-0_0-C/OUTPUT.jld2 $PPFILE BBALANCE.jld2 $RAM/0_14-0_1-0_03-0_0-C &
wait

rmdir $RAM
