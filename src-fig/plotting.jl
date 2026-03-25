# plotting.jl
# Various helper functions for making figures
# -------------------------------------------------------------

# -------------------------------------------------------------
using JLD2
using ImageFiltering: imfilter, Kernel.gaussian
using OffsetArrays: no_offset_view
using Printf
# -------------------------------------------------------------

# -------------------------------------------------------------
const b_step = 0.00004
const b_levels = -0.002:0.00004:0.0004
# -------------------------------------------------------------

# -------------------------------------------------------------
function iterations_times(file)
    iterations = keys(file["timeseries/t"])
    times = [file["timeseries/t/$i"] for i in iterations]

    return iterations, times
end

iterations_times(filename::String) = jldopen(iterations_times, filename)

@inline prettytime(t) = @sprintf "%06.1f" t

simulation_parameters(file) = file["metadata/parameters"]
simulation_parameters(filename::String) = jldopen(simulation_parameters, filename)

@inline function variable_strain_rate(t, sp)
    turnon = max(1-exp(-sp.f * t / 15), 0)
    return sp.α * turnon
end

@inline function surface_b_flux(t, sp) 
    turnon = 1 - exp(-sp.f*(t - sp.start_time) / 20)
    return sp.B * turnon
end

function bin_counts(a, bins) 
    map(bins[1:end-1], bins[2:end]) do bl, br
        sum(bl .< a .< br)
    end
end

function normalized_time(file)
    iterations, times = iterations_times(file)
    sp = simulation_parameters(file)

    ds = map(enumerate(times)) do (i, t)
        Δt = times[i] - times[max(1, i-1)]
        variable_strain_rate(t, sp) * Δt
    end
    
    return cumsum(ds)
end

@inline function interpolate(x, xs, fs)
    i1 = findlast(xs .<= x)
    i1 == nothing && return fs[end]
    i2 = i1 + 1
    i2 > length(fs) && return fs[i1]
    df = (x - xs[i1]) * (fs[i2] - fs[i1]) / (xs[i2] - xs[i1])
    return fs[i1] + df
end
# -------------------------------------------------------------

# -------------------------------------------------------------
include("filt.jl")
include("subfig_label.jl")
include("grid_nodes.jl")
include("bounds.jl")
include("get_field.jl")
include("timeseries_of.jl")
include("time_average_of.jl")
# -------------------------------------------------------------
