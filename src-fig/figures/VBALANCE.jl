VBALANCE_term_labels = (;
    advection_x = L"\text{Across-front advection}",
    advection_background = L"\text{Background advection}",
    advection_z = L"\text{Vertical advection}",
    mixing_x = L"\text{Across-front mixing}",
    mixing_z = L"\text{Vertical mixing}",
    coriolis_y = L"\text{Coriolis}",
    strain_y = L"\text{Strain}",
    sponge = L"\text{Sponge}"
)

function check_VBALANCE(run_id; N_window=1)
    foldername = joinpath(scratchpath, run_id)
    
    suffix = N_window == 1 ? "" : "-$N_window"
    MEAN = joinpath(foldername, "MEAN$(suffix).jld2")
    VBALANCE = joinpath(foldername, "VBALANCE$(suffix).jld2")
        
    iterations, times = iterations_times(VBALANCE)
    sp = simulation_parameters(VBALANCE)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(VBALANCE)
    
    fig = Figure(; size=(1000, 300), fontsize=18)
    ax_actual = Axis(fig[1, 1]; 
        xlabel = x_label,
        ylabel = z_label,
        limits = (-sp.Lh / 2x_unit, sp.Lh / 2x_unit, -sp.Lz / z_unit, 0),
        xticks = [-1, 0, 1]
    )
    
    ax_total = Axis(fig[1 ,2]; 
        xlabel = x_label,
        ylabel = z_label,
        limits = (-sp.Lh / 2x_unit, sp.Lh / 2x_unit, -sp.Lz / z_unit, 0),
        xticks = [-1, 0, 1]
    )
    hideydecorations!(ax_total; ticks=false)
    
    ax_difference = Axis(fig[1 ,3]; 
        xlabel = x_label,
        ylabel = z_label,
        limits = (-sp.Lh / 2x_unit, sp.Lh / 2x_unit, -sp.Lz / z_unit, 0),
        xticks = [-1, 0, 1]
    )
    hideydecorations!(ax_difference; ticks=false)
    
    dvdt_actual = (get_field(MEAN, "v_bar", iterations[end]) .- get_field(MEAN, "v_bar", iterations[end-1])) ./ (times[end] - times[end-1]) ./ (sp.f^2 * sp.L)
    dvdt_total = get_field(VBALANCE, "total", iterations[end]) ./ (sp.f^2 * sp.L)
    dvdt_difference = dvdt_total .- dvdt_actual
    
    colorrange = (-1, 1)
    heatmap!(ax_actual, xsᶜ ./ x_unit, zsᶜ, dvdt_actual; colormap=:balance, colorrange)
    heatmap!(ax_total, xsᶜ ./ x_unit, zsᶜ, dvdt_total; colormap=:balance, colorrange)
    ht = heatmap!(ax_difference, xsᶜ ./ x_unit, zsᶜ, dvdt_difference; colormap=:balance, colorrange)
    
    Colorbar(fig[2, 1:3], ht; vertical=false, flipaxis=false, label=L"A / f^2L_D")
    fig
end

function terms_VBALANCE(run_id, frames, filename;
    fig_kw = NamedTuple(),
    ax_kw = NamedTuple(),
    record_kw = NamedTuple(),
    N_window = 1
    )
    foldername = joinpath(scratchpath, run_id)
    
    suffix = N_window == 1 ? "" : "-$N_window"
    MEAN = joinpath(foldername, "MEAN$(suffix).jld2")
    VBALANCE = joinpath(foldername, "VBALANCE$(suffix).jld2")
    
    sp = simulation_parameters(MEAN)
    iterations, times = iterations_times(MEAN)
    
    n = Observable(frames[1])
    t = @lift interp_time($n, times)
    
    fts_b_bar = FieldTimeSeries(MEAN, "b_bar")
    
    fieldnames = (;
        advection_x = "advection_x",
        advection_background = "advection_background",
        advection_z = "advection_z",
        mixing_x = "mixing_x",
        mixing_z = "mixing_z",
        coriolis_y = "coriolis_y",
        strain_y = "strain_y",
        sponge = "sponge"
    )

    fts = NamedTuple(k => FieldTimeSeries(VBALANCE, v; backend=OnDisk()) for (k, v) in pairs(fieldnames))
    
    title = @lift let t_hr = @sprintf "%.0f" ($t / 3600)
        L"\text{Terms in }\overline{v}\text{ balance}\quad t = %$t_hr \, \text{hr}"
    end
    
    b_bar = @lift nov(fts_b_bar[Time($t)][:, 1, :]) ./ sp.Δb

    field_data = NamedTuple(k => @lift nov(v[Time($t)][:, 1, :]) ./ (sp.L * sp.f^2) for (k, v) in pairs(fts))
    
    fig = Figure(; 
        size=(1000, 600),
        fig_kw...
    )
    Label(fig[1, 1:4], title)
    
    ax_kw = (;
        xlabel = x_label,
        ylabel = z_label,
        limits = (-sp.Lh / 2x_unit, sp.Lh / 2x_unit, -sp.Lz / z_unit, 0),
        xticks = [-1, 0, 1]
    )

    axes = (;
        advection_x = Axis(fig[2, 1]; ax_kw..., title=VBALANCE_term_labels.advection_x),
        advection_background = Axis(fig[2, 2]; ax_kw..., title=VBALANCE_term_labels.advection_background),
        advection_z = Axis(fig[2, 3]; ax_kw..., title=VBALANCE_term_labels.advection_z),
        mixing_x = Axis(fig[3, 1]; ax_kw..., title=VBALANCE_term_labels.mixing_x),
        mixing_z = Axis(fig[3, 2]; ax_kw..., title=VBALANCE_term_labels.mixing_z),
        coriolis_y = Axis(fig[2, 4]; ax_kw..., title=VBALANCE_term_labels.coriolis_y),
        strain_y = Axis(fig[3, 3]; ax_kw..., title=VBALANCE_term_labels.strain_y),
        sponge = Axis(fig[3, 4]; ax_kw..., title=VBALANCE_term_labels.sponge),
    )
    

    hidexdecorations!(axes.advection_x; ticks=false)
    hidexdecorations!(axes.advection_background; ticks=false)
    hidexdecorations!(axes.advection_z; ticks=false)
    hidexdecorations!(axes.coriolis_y; ticks=false)

    hideydecorations!(axes.advection_background; ticks=false)
    hideydecorations!(axes.advection_z; ticks=false)
    hideydecorations!(axes.coriolis_y; ticks=false)
    hideydecorations!(axes.mixing_z; ticks=false)
    hideydecorations!(axes.strain_y; ticks=false)
    hideydecorations!(axes.sponge; ticks=false)

    hts = NamedTuple(
        map(keys(axes)) do k
            xs = nov(xnodes(fts[k]; with_halos=true)) ./ x_unit
            zs = nov(znodes(fts[k]; with_halos=true)) ./ z_unit
            data = field_data[k]
            colormap = :balance
            colorrange = (-1, 1)
    
            k => heatmap!(axes[k], xs, zs, data; colormap, colorrange)
        end
    )
    
    map(keys(axes)) do k
        xs = nov(xnodes(fts_b_bar; with_halos=true)) ./ x_unit
        zs = nov(znodes(fts_b_bar; with_halos=true)) ./ z_unit
        data = b_bar
        levels = b_levels(fts_b_bar, sp) ./ sp.Δb
        color = (:black, 0.5)

        contour!(axes[k], xs, zs, data; levels, color)
    end

    Colorbar(fig[2:3, 5], hts.advection_x; label=L"\text{Term} / f^2L_D")

    prettyrecord(n, fig, filename, frames; record_kw...)

    return fig
end