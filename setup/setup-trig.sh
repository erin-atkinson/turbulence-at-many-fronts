#!/bin/bash
# Run this in a login node to install Julia for Trillium (GPU)

module load julia/1.12.5
export JULIA_DEPOT_PATH=$SCRATCH/julia-trig

julia -t 4 -- setup/setup-trig.jl
