#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=0:45:00
#SBATCH --job-name=cooling-only-SLICES
#SBATCH --output=../scratch/logs/cooling-only-SLICES.txt

module load julia/1.12.5

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

PPFILE=SLICES
RAM=/dev/shm/cooling-only-SLICES
mkdir $RAM
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_01-0_0-C $PPFILE $RAM/0_1-0_1-0_01-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_02-0_0-C $PPFILE $RAM/0_1-0_1-0_02-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_03-0_0-C $PPFILE $RAM/0_1-0_1-0_03-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_04-0_0-C $PPFILE $RAM/0_1-0_1-0_04-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_05-0_0-C $PPFILE $RAM/0_1-0_1-0_05-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_06-0_0-C $PPFILE $RAM/0_1-0_1-0_06-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_07-0_0-C $PPFILE $RAM/0_1-0_1-0_07-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_08-0_0-C $PPFILE $RAM/0_1-0_1-0_08-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_09-0_0-C $PPFILE $RAM/0_1-0_1-0_09-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_1-0_0-C $PPFILE $RAM/0_1-0_1-0_1-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_11-0_0-C $PPFILE $RAM/0_1-0_1-0_11-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_12-0_0-C $PPFILE $RAM/0_1-0_1-0_12-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_13-0_0-C $PPFILE $RAM/0_1-0_1-0_13-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_14-0_0-C $PPFILE $RAM/0_1-0_1-0_14-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_15-0_0-C $PPFILE $RAM/0_1-0_1-0_15-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_2-0_0-C $PPFILE $RAM/0_1-0_1-0_2-0_0-C &
wait

rm $RAM -rf
