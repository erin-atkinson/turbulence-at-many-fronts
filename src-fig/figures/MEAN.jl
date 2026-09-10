@doc raw"""
    mean_fields(run_id, frames, filename=nothing;
        record_kw = NamedTuple(),
        N_window = 1
    )
Create and animate a figure of along-front velocity, buoyancy and streamfunction
"""
function mean_fields(run_id, frames, filename=nothing;
        record_kw = NamedTuple(),
        N_window = 1
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
        xs = nov(xnodes(fts_b_bar; with_halos=true)) ./ 1000
        zs = nov(znodes(fts_b_bar; with_halos=true))
        data = b_bar
        levels = b_levels(fts_b_bar, sp) ./ sp.Δb
        color = (:black, 0.5)

        contour!(ax_u, xs, zs, data; levels, color)
        contour!(ax_v, xs, zs, data; levels, color)
        contour!(ax_w, xs, zs, data; levels, color)
    end

    Colorbar(fig[3, 1], ht_u; flipaxis=false, vertical=false, label=background ? tot_u_bar_label : u_bar_label)
    Colorbar(fig[3, 2], ht_v; flipaxis=false, vertical=false, label=v_bar_label)
    Colorbar(fig[3, 3], ht_w; flipaxis=false, vertical=false, label=w_bar_label)

    colgap!(fig.layout, 40)
    prettyrecord(n, fig, filename, frames; record_kw...)

    return fig
end

@doc raw"""
    mean_fields(run_id, filename;
        record_kw = NamedTuple(),
        N_window = 1
    )
Create a hovmoller plot of total across-front velocity, along-front velocity and vertical velocity 
"""
function mean_hovmoller(foldername, z;
    fig_kw = NamedTuple(),
    ax_kw = NamedTuple(),
    background = true,
    remove_mean = false,
    N_window = 1
    )

    suffix = N_window == 1 ? "" : "-$N_window"
    MEAN = joinpath(foldername, "MEAN$suffix.jld2")

    fts_u_bar = FieldTimeSeries(MEAN, "u_bar")
    fts_v_bar = FieldTimeSeries(MEAN, "v_bar")
    fts_w_bar = FieldTimeSeries(MEAN, "w_bar")
    fts_b_bar = FieldTimeSeries(MEAN, "b_bar")
    
    sp = simulation_parameters(MEAN)
    times = fts_u_bar.times
    
    n = Observable(frames[1])
    t = @lift interp_time($n, times)
    
    title = let z_str = @sprintf "%.0f" z
        L"\text{Mean fields}\quad z = %$z_stw \, \text{m}"
    end
    
    U = @lift if background
        [velocity_profile(x, sp) * variable_strain_rate($t, sp) for x in xnodes(fts_u_bar; with_halos=true), z in 1:1]
    else
        0
    end
    
    u_bar = @lift nov(fts_u_bar[Time($t)][:, 1, :] .+ $U) .* 100 
    v_bar = @lift nov(fts_v_bar[Time($t)][:, 1, :]) .* 100
    w_bar = @lift nov(fts_w_bar[Time($t)][:, 1, :]) .* 1000
    b_bar = @lift nov(fts_b_bar[Time($t)][:, 1, :]) ./ sp.Δb
    
    fig = Figure(; 
        size=(1000, 400),
        fig_kw...
    )
    Label(fig[1, 1:3], title)
    
    ax_kw = (;
        xlabel = L"x / \text{km}",
        ylabel = L"z / \text{m}",
        limits = (-sp.Lh / 2000, sp.Lh / 2000, -sp.Lz, 0)
    )

    ax_u = Axis(fig[2, 1]; ax_kw...)
    ax_v = Axis(fig[2, 2]; ax_kw...)
    ax_w = Axis(fig[2, 3]; ax_kw...)

    hideydecorations!(ax_v; ticks=false)
    hideydecorations!(ax_w; ticks=false)

    ht_u = begin
        xs = nov(xnodes(fts_u_bar; with_halos=true)) ./ 1000
        zs = nov(znodes(fts_u_bar; with_halos=true))
        data = u_bar
        colormap = :balance
        colorrange = (-10, 10)

        heatmap!(ax_u, xs, zs, data; colormap, colorrange)
    end

    ht_v = begin
        xs = nov(xnodes(fts_v_bar; with_halos=true)) ./ 1000
        zs = nov(znodes(fts_v_bar; with_halos=true))
        data = v_bar
        colormap = :balance
        colorrange = (-10, 10)

        heatmap!(ax_v, xs, zs, data; colormap, colorrange)
    end

    ht_w = begin
        xs = nov(xnodes(fts_w_bar; with_halos=true)) ./ 1000
        zs = nov(znodes(fts_w_bar; with_halos=true))
        data = w_bar
        colormap = :balance
        colorrange = (-10, 10)

        heatmap!(ax_w, xs, zs, data; colormap, colorrange)
    end

    begin 
        xs = nov(xnodes(fts_b_bar; with_halos=true)) ./ 1000
        zs = nov(znodes(fts_b_bar; with_halos=true))
        data = b_bar
        levels = b_levels(fts_b_bar, sp) ./ sp.Δb
        color = (:black, 0.5)

        contour!(ax_u, xs, zs, data; levels, color)
        contour!(ax_v, xs, zs, data; levels, color)
        contour!(ax_w, xs, zs, data; levels, color)
    end

    Colorbar(fig[3, 1], ht_u; flipaxis=false, vertical=false, label=background ? tot_u_bar_label : u_bar_label)
    Colorbar(fig[3, 2], ht_v; flipaxis=false, vertical=false, label=v_bar_label)
    Colorbar(fig[3, 3], ht_w; flipaxis=false, vertical=false, label=w_bar_label)

    colgap!(fig.layout, 40)
    prettyrecord(n, fig, filename, frames; record_kw...)

    return fig
end
