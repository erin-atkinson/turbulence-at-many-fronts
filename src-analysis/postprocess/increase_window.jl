using Oceananigans
using JLD2
using Printf

using Oceananigans.Fields: compute_at!
using Oceananigans.OutputWriters: saveproperty!
using Oceananigans: fill_halo_regions!

include("update_clock.jl")
include("update_fields.jl")

t0 = time()

# Simulation output file
RAW = ARGS[1]
foldername = dirname(RAW)
filename = basename(RAW)

# Number of iterations to sum 
N_window = parse(Int64, ARGS[2])
outputname = splitext(filename)[1] * "-$N_window"

# Possible third argument is a temporary location
buffer = length(ARGS) > 2 ? ARGS[3] : foldername
mkpath(buffer)

# Path to output
BUFFER = joinpath(buffer, "$outputname.jld2")

# Final output folder
PROCESSED = joinpath(foldername, "$outputname.jld2")

# Read simulation state
fds = FieldDataset(RAW; backend=OnDisk())

if "parameters" ∉ keys(fds.metadata)
    parameterfilename = joinpath(foldername, "parameters.jld2")
    fds.metadata["parameters"] = jldopen(file->file["simulation"], parameterfilename)
end

fieldnames = Symbol.(keys(fds.fields))
const sp = fds.metadata["parameters"]

# Setup grid, times and iterations
grid = fds.u.grid
times = fds.u.times
iterations = jldopen(file->keys(file["timeseries/t"]), RAW)
iterations = parse.(Int, iterations)
frames = 1:N_window:length(iterations)

# Named tuple of current simulation state fields
rawfields = NamedTuple(k => deepcopy(fds[k][1]) for k in fieldnames)
input_fields = rawfields

# Output fields
output_fields = NamedTuple(k => deepcopy(v) for (k, v) in pairs(rawfields))

# Accumulation fields
accfields = NamedTuple(k => Field(output_fields[k] + rawfields[k] / N_window) for k in fieldnames)

output_fds = FieldDataset(times, output_fields; 
    backend = OnDisk(), 
    path = BUFFER,
    metadata = fds.metadata
)

@info "Computing..."
t1 = time()
t2 = t1
for (i, frame) in enumerate(frames)
    iteration = iterations[frame]
    t = times[frame]
    
    for k in fieldnames
        set!(rawfields[k], 0)
        set!(output_fields[k], 0)
        set!(accfields[k], 0)
    end

    for n in 0:(N_window - 1)
        map(k->set!(rawfields[k], fds[k][frame + n]), fieldnames)
        map(k->compute_at!(accfields[k], frame + n), fieldnames)
        map(k->set!(output_fields[k], accfields[k]), fieldnames)
    end
    
    set!(output_fds, iteration, t; output_fields...)

    # Little bit of timekeeping for convenience
    
    tstr = if i < 2
        global t2 = time()
        setup_time = t2 - t0
        tstr = @sprintf "Setup: %.2f s" setup_time
    else
        t3 = time()
        setup_time = t2 - t0
        avg_time = (t3 - t2) / (i - 10)
        total_time = setup_time + avg_time * (length(frames) - 10)
        @sprintf "Setup: %.2f s, Avg: %.2f s, Total: %.2f s" setup_time avg_time total_time
    end
    
    print("$(frames[1]) -> $frame -> $(frames[end]) | $tstr\r")
end
println()

# Write grid to file
jldopen(file->saveproperty!(file, "grid", grid), BUFFER, "a")

if !isequal(BUFFER, PROCESSED)
    @info "Moving from $BUFFER to $PROCESSED"
    mv(BUFFER, PROCESSED; force=true)
end

@info "Done!"