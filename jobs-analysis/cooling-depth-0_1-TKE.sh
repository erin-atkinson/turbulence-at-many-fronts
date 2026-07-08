#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=2:30:00
#SBATCH --job-name=cooling-depth-0_1-TKE
#SBATCH --output=../scratch/logs/cooling-depth-0_1-TKE.txt

module load julia/1.12.5

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

PPFILE=TKE
RAM=/dev/shm/cooling-depth-0_1-TKE
mkdir $RAM
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_06-0_1-0_01-0_0-C $PPFILE $RAM/0_06-0_1-0_01-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_06-0_1-0_03-0_0-C $PPFILE $RAM/0_06-0_1-0_03-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_06-0_1-0_05-0_0-C $PPFILE $RAM/0_06-0_1-0_05-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_06-0_1-0_1-0_0-C $PPFILE $RAM/0_06-0_1-0_1-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_06-0_1-0_2-0_0-C $PPFILE $RAM/0_06-0_1-0_2-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_01-0_0-C $PPFILE $RAM/0_1-0_1-0_01-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_03-0_0-C $PPFILE $RAM/0_1-0_1-0_03-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_05-0_0-C $PPFILE $RAM/0_1-0_1-0_05-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_1-0_0-C $PPFILE $RAM/0_1-0_1-0_1-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-0_1-0_2-0_0-C $PPFILE $RAM/0_1-0_1-0_2-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_14-0_1-0_01-0_0-C $PPFILE $RAM/0_14-0_1-0_01-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_14-0_1-0_03-0_0-C $PPFILE $RAM/0_14-0_1-0_03-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_14-0_1-0_05-0_0-C $PPFILE $RAM/0_14-0_1-0_05-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_14-0_1-0_1-0_0-C $PPFILE $RAM/0_14-0_1-0_1-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_14-0_1-0_2-0_0-C $PPFILE $RAM/0_14-0_1-0_2-0_0-C &
wait

rm $RAM -rf
