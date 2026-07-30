# plotting.jl
# Various helper functions for making figures
# -------------------------------------------------------------

# -------------------------------------------------------------
using Oceananigans
using CairoMakie
using CairoMakie.Colors
using Statistics
using Printf
using JLD2
using ImageFiltering: imfilter, Kernel.gaussian
using OffsetArrays: no_offset_view
using Printf
using Oceananigans.Units: Time
# -------------------------------------------------------------

# -------------------------------------------------------------
# Constants used for plotting 
scratchpath = "/home/atkin163/scratch/turbulence-at-many-fronts"

# Labels
const u_bar_label = L"\overline{u} / \text{cm}\,\text{s}^{-1}"
const tot_u_bar_label = L"(\overline{u} + U) / \text{cm}\,\text{s}^{-1}"
const v_bar_label = L"\overline{v} / \text{cm}\,\text{s}^{-1}"
const w_bar_label = L"\overline{w} / \text{mm}\,\text{s}^{-1}"
const b_bar_label = L"\overline{b} / \text{m}\,\text{s}^{-2}" 

const x_label = L"x / \text{km}"
const y_label = L"y / \text{km}"
const z_label = L"z / \text{m}"
const t_label = L"t / \text{hr}"

const u_unit = 1e-2
const v_unit = 1e-2
const w_unit = 1e-3

const x_unit = 1e3
const y_unit = 1e3
const z_unit = 1
const t_unit = 3600

# Distance between b contours relative to Δb
const b_step = 1/6
b_levels(fds::FieldDataset) = minimum(interior(fds.b_bar[end], :, 1, :)):(fds.metadata["parameters"].Δb * b_step):maximum(interior(fds.b_bar[1], :, 1, :))
b_levels(fts, sp) = minimum(interior(fts[end], :, 1, :)):(sp.Δb * b_step):maximum(interior(fts[1], :, 1, :))
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

@inline function interp_time(n, times)
    i = Integer(floor(n))
    ξ = n - i
    return times[i] + ξ * (times[min(i+1, length(times))] - times[i])
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
include("front_width.jl")
include("record.jl")
include("filenames.jl")
include("../src-analysis/terms/forcing_bc_funcs.jl")
# -------------------------------------------------------------
