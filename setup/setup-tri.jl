using Pkg
@info "Setting up Julia depot for Trillium at $DEPOT_PATH"

Pkg.add(; name="Oceananigans", version="0.107.4")
Pkg.pin("Oceananigans")

Pkg.add(; name="JLD2", version="0.6.4")
Pkg.pin("JLD2")

Pkg.add("OffsetArrays")
Pkg.add("SpecialFunctions")

try
    Pkg.precompile()
catch e
    println("First precompile failed, retrying")
    Pkg.precompile()
end

Pkg.status()
