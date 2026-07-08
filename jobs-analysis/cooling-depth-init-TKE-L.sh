#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=0:45:00
#SBATCH --job-name=cooling-depth-init-TKE-L
#SBATCH --output=../scratch/logs/cooling-depth-init-TKE-L.txt

module load julia/1.12.5

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

PPFILE=TKE-L
RAM=/dev/shm/cooling-depth-init-TKE-L
mkdir $RAM
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_06-init-0_01-0_0-C $PPFILE $RAM/0_06-init-0_01-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_06-init-0_03-0_0-C $PPFILE $RAM/0_06-init-0_03-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_06-init-0_05-0_0-C $PPFILE $RAM/0_06-init-0_05-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_06-init-0_1-0_0-C $PPFILE $RAM/0_06-init-0_1-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_06-init-0_2-0_0-C $PPFILE $RAM/0_06-init-0_2-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-init-0_01-0_0-C $PPFILE $RAM/0_1-init-0_01-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-init-0_03-0_0-C $PPFILE $RAM/0_1-init-0_03-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-init-0_05-0_0-C $PPFILE $RAM/0_1-init-0_05-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-init-0_1-0_0-C $PPFILE $RAM/0_1-init-0_1-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_1-init-0_2-0_0-C $PPFILE $RAM/0_1-init-0_2-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_14-init-0_01-0_0-C $PPFILE $RAM/0_14-init-0_01-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_14-init-0_03-0_0-C $PPFILE $RAM/0_14-init-0_03-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_14-init-0_05-0_0-C $PPFILE $RAM/0_14-init-0_05-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_14-init-0_1-0_0-C $PPFILE $RAM/0_14-init-0_1-0_0-C &
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/0_14-init-0_2-0_0-C $PPFILE $RAM/0_14-init-0_2-0_0-C &
wait

rm $RAM -rf
