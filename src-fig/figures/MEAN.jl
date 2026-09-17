@doc raw"""
    mean_figure(run_id, frames;
        record_kw = NamedTuple(),
        N_window = 1
        filename = joinpath(run_id, "$run_id-mean")
    )
Create and animate a figure of along-front velocity, buoyancy and streamfunction
"""
function mean_figure(run_id, frames;
        record_kw = NamedTuple(),
        N_window = 1,
        filename = joinpath(run_id, "$run_id-mean")
    )
    
    MEAN = filepath(run_id, "MEAN", N_window)

    fts_v_bar = FieldTimeSeries(MEAN, "v_bar")
    fts_b_bar = FieldTimeSeries(MEAN, "b_bar")
    fts_ψ = FieldTimeSeries(MEAN, "ψ")
    
    sp = simulation_parameters(MEAN)
    times = fts_v_bar.times
    
    n = Observable(frames[1])
    t = @lift interp_time($n, times)
    
    title = @lift let t_hr = @sprintf "%.0f" ($t / 3600)
        L"\text{Mean fields}\quad t = %$t_hr \, \text{hr}"
    end
    
    v_bar = @lift nov(fts_v_bar[Time($t)][:, 1, :]) ./ v_unit
    b_bar = @lift nov(fts_b_bar[Time($t)][:, 1, :]) ./ sp.Δb
    ψ = @lift nov(fts_ψ[Time($t)][:, 1, :]) ./ ψ_unit
    
    fig = Figure(; size=(figure_width, 400), fontsize)
    Label(fig[1, 1:3], title)
    
    ax_kw = (;
        xlabel = x_label,
        ylabel = z_label,
        limits = transect_limits(sp)
    )

    ax_v = Axis(fig[2, 1]; ax_kw...)
    ax_b = Axis(fig[2, 2]; ax_kw...)
    ax_ψ = Axis(fig[2, 3]; ax_kw...)

    hideydecorations!(ax_b; ticks=false)
    hideydecorations!(ax_ψ; ticks=false)

    ht_v = begin
        xs = nov(xnodes(fts_v_bar; with_halos=true)) ./ x_unit
        zs = nov(znodes(fts_v_bar; with_halos=true))
        data = v_bar
        colormap = :balance
        colorrange = (-10, 10)

        heatmap!(ax_v, xs, zs, data; colormap, colorrange)
    end

    ht_b = begin
        xs = nov(xnodes(fts_v_bar; with_halos=true)) ./ x_unit
        zs = nov(znodes(fts_v_bar; with_halos=true))
        data = b_bar
        colormap = :balance
        colorrange = (-10, 10)

        heatmap!(ax_b, xs, zs, data; colormap, colorrange)
    end

    ht_ψ = begin
        xs = nov(xnodes(fts_ψ_bar; with_halos=true)) ./ x_unit
        zs = nov(znodes(fts_ψ_bar; with_halos=true))
        data = ψ
        colormap = :balance
        colorrange = (-10, 10)

        heatmap!(ax_ψ, xs, zs, data; colormap, colorrange)
    end

    begin 
        xs = nov(xnodes(fts_b_bar; with_halos=true)) ./ x_unit
        zs = nov(znodes(fts_b_bar; with_halos=true))
        data = b_bar
        levels = b_levels(fts_b_bar, sp) ./ sp.Δb
        color = (:black, 0.5)

        contour!(ax_u, xs, zs, data; levels, color)
        contour!(ax_v, xs, zs, data; levels, color)
        contour!(ax_w, xs, zs, data; levels, color)
    end

    Colorbar(fig[3, 1], ht_v; flipaxis=false, vertical=false, label=v_bar_label)
    Colorbar(fig[3, 2], ht_b; flipaxis=false, vertical=false, label=b_bar_label)
    Colorbar(fig[3, 3], ht_ψ; flipaxis=false, vertical=false, label=ψ_label)

    colgap!(fig.layout, 40)
    prettyrecord(n, fig, filename, frames; record_kw...)

    return fig
end

function build_hovmoller(filename, field, frames, z)
    println("Slicing $field from $filename at z=$z")
    fts = FieldTimeSeries(filename, field)

    times = fts.u_bar.times
    t = interp_time(frames[1], times)

    c = similar(fts[Time(t)])
    c_slice = Field(ZSlice(c, z))

    c_hovmoller = zeros(eltype(c_slice), length(frames), size(c_slice, 1))

    for (n, frame) in enumerate(frames)
        t = interp_time(frame, times)
        set!(c, fts[Time(t)])

        compute!(c_slice)
        fill_halo_regions!(c_slice)

        c_hovmoller[:, n] .= interior(c_slice, :, 1, 1)
    end

    return c_hovmoller
end

@doc raw"""
    mean_hovmoller(run_id, frames, z;
        N_window = 1
        filename = joinpath(imagepath, run_id, "$run_id-mean_hovmoller.png"),
        background = true
    )
Create hovmoller plots of total across-front velocity, along-front velocity and vertical velocity with buoyancy contours
"""
function mean_hovmoller(run_id, frames, z;
    fig_kw = NamedTuple(),
    ax_kw = NamedTuple(),
    background = true,
    N_window = 1,
    filename = joinpath(imagepath, run_id, "$run_id-mean_hovmoller.png")
    )

    MEAN = filepath(run_id, "MEAN", N_window)

    sp = simulation_parameters(MEAN)
    iterations, times = iteration_times(MEAN)
    times = [interp_time(frame, times) for frame in frames]

    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(MEAN)

    U = if background
        [velocity_profile(x, sp) * variable_strain_rate(t, sp) for x in xsᶠ, t in times]
    else
        0
    end

    u_hovmoller = build_hovmoller(MEAN, "u_bar", frames, z) .+ U
    v_hovmoller = build_hovmoller(MEAN, "v_bar", frames, z)
    w_hovmoller = build_hovmoller(MEAN, "w_bar", frames, z)
    b_hovmoller = build_hovmoller(MEAN, "b_bar", frames, z)

    title = let z_str = @sprintf "%.0f" z
        L"\text{Mean fields}\quad z = %$z_str \, \text{m}"
    end
    
    fig = Figure(; size=(figure_width, 600), fontsize, fig_kw...)
    Label(fig[1, 1:3], title)
    
    ax_kw = (;
        xlabel = x_label,
        ylabel = t_label,
        limits = (-sp.Lh / 2x_unit, sp.Lh / 2x_unit, times[1] / t_unit, times[end] / t_unit),
        ax_kw...
    )
    
    ax_u = Axis(fig[2, 1]; ax_kw...)
    ax_v = Axis(fig[2, 2]; ax_kw...)
    ax_w = Axis(fig[2, 3]; ax_kw...)

    hideydecorations!(ax_v; ticks=false)
    hideydecorations!(ax_w; ticks=false)

    ht_u = begin
        xs = xsᶠ ./ x_unit
        zs = times ./ t_unit
        data = u_hovmoller ./ u_unit
        colormap = :balance
        colorrange = (-10, 10)

        heatmap!(ax_u, xs, zs, data; colormap, colorrange)
    end

    ht_v = begin
        xs = xsᶜ ./ x_unit
        zs = times ./ t_unit
        data = v_hovmoller ./ v_unit
        colormap = :balance
        colorrange = (-10, 10)

        heatmap!(ax_v, xs, times, data; colormap, colorrange)
    end

    ht_w = begin
        xs = xsᶜ ./ x_unit
        zs = times ./ t_unit
        data = w_hovmoller ./ w_unit
        colormap = :balance
        colorrange = (-10, 10)

        heatmap!(ax_w, xs, zs, data; colormap, colorrange)
    end

    begin 
        xs = xsᶜ ./ x_unit
        zs = times ./ t_unit
        data = b_hovmoller ./ sp.Δb

        levels = minimum(data):(1/6):maximum(data)
        color = (:black, 0.5)

        contour!(ax_u, xs, zs, data; levels, color)
        contour!(ax_v, xs, zs, data; levels, color)
        contour!(ax_w, xs, zs, data; levels, color)
    end

    Colorbar(fig[3, 1], ht_u; flipaxis=false, vertical=false, label=background ? tot_u_bar_label : u_bar_label)
    Colorbar(fig[3, 2], ht_v; flipaxis=false, vertical=false, label=v_bar_label)
    Colorbar(fig[3, 3], ht_w; flipaxis=false, vertical=false, label=w_bar_label)

    colgap!(fig.layout, 40)
    save(filename, fig; record_kw...)

    return fig
end
