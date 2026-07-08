#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=3:00:00
#SBATCH --job-name=ppBBALANCE
#SBATCH --output=../scratch/logs/ppBBALANCE.txt

module load julia/1.10.10

# Launch from scratch
export JULIA_DEPOT_PATH=$SCRATCH/julia-tri
cd ~/turbulence-at-many-fronts

PPFILE=BBALANCE
RAM=/dev/shm/$PPFILE
mkdir $RAM

a=0_30
b=1_0
c=1_0
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/a$a-b$b-C$c $PPFILE $RAM/a$a-b$b-C$c &

c=1_4
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/a$a-b$b-C$c $PPFILE $RAM/a$a-b$b-C$c &

c=1_8
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/a$a-b$b-C$c $PPFILE $RAM/a$a-b$b-C$c &

c=2_0
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/a$a-b$b-C$c $PPFILE $RAM/a$a-b$b-C$c &

#c=2_2
#julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/a$a-b$b-C$c $PPFILE $RAM/a$a-b$b-C$c &

c=2_6
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/a$a-b$b-C$c $PPFILE $RAM/a$a-b$b-C$c &

c=3_0
julia -t 24 -- src-analysis/postprocess/postprocess.jl $SCRATCH/turbulence-at-many-fronts/a$a-b$b-C$c $PPFILE $RAM/a$a-b$b-C$c &

wait
rm $RAM -rf