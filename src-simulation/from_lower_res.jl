function need_to_interpolate(output_folder, initialised_folder)
    initialised_folder == "" && return false
    # Check the output folder for checkpoints
    any(x -> occursin(r"^checkpoint", x), readdir(output_folder)) && return false
    any(x -> occursin("OUTPUT.jld2", x), readdir(initialised_folder)) && return true

    return false
end

function set_interpolated_init!(model, initialised_folder, sp)
    OUTPUT = joinpath(initialised_folder, "OUTPUT.jld2")
    
    u_fts = FieldTimeSeries(OUTPUT, "u"; backend=OnDisk())
    v_fts = FieldTimeSeries(OUTPUT, "v"; backend=OnDisk())
    w_fts = FieldTimeSeries(OUTPUT, "w"; backend=OnDisk())
    b_fts = FieldTimeSeries(OUTPUT, "b"; backend=OnDisk())

    u_lower = Field(Average(u_fts[Oceananigans.Units.Time(sp.start_time)]; dims=2))
    v_lower = Field(Average(v_fts[Oceananigans.Units.Time(sp.start_time)]; dims=2))
    w_lower = Field(Average(w_fts[Oceananigans.Units.Time(sp.start_time)]; dims=2))
    b_lower = Field(Average(b_fts[Oceananigans.Units.Time(sp.start_time)]; dims=2))

    # Only have to do this once so maybe this is okay...
    @info "Interpolating u"
    set!(model, u=(x, y, z) -> Oceananigans.interpolate((x, z), u_lower))
    @info "Interpolating v"
    set!(model, v=(x, y, z) -> Oceananigans.interpolate((x, z), v_lower))
    @info "Interpolating w"
    set!(model, w=(x, y, z) -> Oceananigans.interpolate((x, z), w_lower))
    @info "Interpolating b"
    set!(model, b=(x, y, z) -> Oceananigans.interpolate((x, z), b_lower))
    
    return nothing
end
