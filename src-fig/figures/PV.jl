# PV.jl

@doc raw"""
    potential_vorticity_figure(run_id, frames;
        record_kw = NamedTuple(),
        N_window = 1,
        filename = joinpath(run_id, "$run_id-pv")
    )
Create a figure of the potential vorticity and its transport
"""
function potential_vorticity_figure(run_id, frames;
    record_kw = NamedTuple(),
    N_window = 1,
    filename = joinpath(run_id, "$run_id-pv")
    )
    
    MEAN = filepath(run_id, "MEAN", N_window)

    fts_v_bar = FieldTimeSeries(MEAN, "v_bar")
    fts_b_bar = FieldTimeSeries(MEAN, "b_bar")
    fts_ψ = FieldTimeSeries(MEAN, "ψ")
    
end
@doc raw"""
    potential_vorticity_vertical_transport(run_id, frames, filename=nothing;
        record_kw = NamedTuple(),
        N_window = 1
    )
Create a 2D timeseries of potential vorticity vertical transport
"""
