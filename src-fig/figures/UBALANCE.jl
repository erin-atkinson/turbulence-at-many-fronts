UBALANCE_term_labels = (;
    advection_x = L"\text{Across-front advection}",
    advection_background = L"\text{Background advection}",
    advection_z = L"\text{Vertical advection}",
    mixing_x = L"\text{Across-front mixing}",
    mixing_z = L"\text{Vertical mixing}",
    coriolis_x = L"\text{Coriolis}",
    strain_x = L"\text{Strain}",
    pressure_x = L"\text{Pressure}",
    sponge = L"\text{Sponge}"
)

function check_UBALANCE(run_id; N_window=1)
    foldername = joinpath(scratchpath, run_id)
    
    suffix = N_window == 1 ? "" : "-$N_window"
    MEAN = joinpath(foldername, "MEAN$(suffix).jld2")
    UBALANCE = joinpath(foldername, "UBALANCE$(suffix).jld2")
        
    iterations, times = iterations_times(UBALANCE)
    sp = simulation_parameters(UBALANCE)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(UBALANCE)
    
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
    
    dudt_actual = (get_field(MEAN, "u_bar", iterations[end]) .- get_field(MEAN, "u_bar", iterations[end-1])) ./ (times[end] - times[end-1]) ./ (sp.f^2 * sp.L)
    dudt_total = get_field(UBALANCE, "total", iterations[end]) ./ (sp.f^2 * sp.L)
    dudt_difference = dudt_total .- dudt_actual
    
    colorrange = (-1, 1)
    heatmap!(ax_actual, xsᶠ ./ x_unit, zsᶜ, dudt_actual; colormap=:balance, colorrange)
    heatmap!(ax_total, xsᶠ ./ x_unit, zsᶜ, dudt_total; colormap=:balance, colorrange)
    ht = heatmap!(ax_difference, xsᶠ ./ x_unit, zsᶜ, dudt_difference; colormap=:balance, colorrange)
    
    Colorbar(fig[2, 1:3], ht; vertical=false, flipaxis=false, label=L"A / f^2L_D")
    fig
end

function terms_UBALANCE(run_id, frames, filename;
    fig_kw = NamedTuple(),
    ax_kw = NamedTuple(),
    record_kw = NamedTuple(),
    N_window = 1
    )
    foldername = joinpath(scratchpath, run_id)
    
    suffix = N_window == 1 ? "" : "-$N_window"
    MEAN = joinpath(foldername, "MEAN$(suffix).jld2")
    UBALANCE = joinpath(foldername, "UBALANCE$(suffix).jld2")
    
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
        coriolis_x = "coriolis_x",
        strain_x = "strain_x",
        pressure_x = "pressure_x",
        sponge = "sponge"
    )

    fts = NamedTuple(k => FieldTimeSeries(UBALANCE, v; backend=OnDisk()) for (k, v) in pairs(fieldnames))
    
    title = @lift let t_hr = @sprintf "%.0f" ($t / 3600)
        L"\text{Terms in }\overline{u}\text{ balance}\quad t = %$t_hr \, \text{hr}"
    end
    
    b_bar = @lift nov(fts_b_bar[Time($t)][:, 1, :]) ./ sp.Δb

    field_data = NamedTuple(k => @lift nov(v[Time($t)][:, 1, :]) ./ (sp.L * sp.f^2) for (k, v) in pairs(fts))
    
    fig = Figure(; 
        size=(1000, 1000),
        fig_kw...
    )
    Label(fig[1, 1:3], title)
    
    ax_kw = (;
        xlabel = x_label,
        ylabel = z_label,
        limits = (-sp.Lh / 2x_unit, sp.Lh / 2x_unit, -sp.Lz / z_unit, 0),
        xticks = [-1, 0, 1]
    )

    axes = (;
        advection_x = Axis(fig[2, 1]; ax_kw..., title=UBALANCE_term_labels.advection_x),
        advection_background = Axis(fig[2, 2]; ax_kw..., title=UBALANCE_term_labels.advection_background),
        advection_z = Axis(fig[2, 3]; ax_kw..., title=UBALANCE_term_labels.advection_z),
        mixing_x = Axis(fig[3, 1]; ax_kw..., title=UBALANCE_term_labels.mixing_x),
        mixing_z = Axis(fig[3, 2]; ax_kw..., title=UBALANCE_term_labels.mixing_z),
        coriolis_x = Axis(fig[4, 1]; ax_kw..., title=UBALANCE_term_labels.coriolis_x),
        pressure_x = Axis(fig[4, 2]; ax_kw..., title=UBALANCE_term_labels.pressure_x),
        strain_x = Axis(fig[3, 3]; ax_kw..., title=UBALANCE_term_labels.strain_x),
        sponge = Axis(fig[4, 3]; ax_kw..., title=UBALANCE_term_labels.sponge),
    )
    

    hidexdecorations!(axes.advection_x; ticks=false)
    hidexdecorations!(axes.advection_background; ticks=false)
    hidexdecorations!(axes.advection_z; ticks=false)
    hidexdecorations!(axes.mixing_x; ticks=false)
    hidexdecorations!(axes.mixing_z; ticks=false)
    hidexdecorations!(axes.strain_x; ticks=false)

    hideydecorations!(axes.advection_background; ticks=false)
    hideydecorations!(axes.advection_z; ticks=false)
    hideydecorations!(axes.mixing_z; ticks=false)
    hideydecorations!(axes.strain_x; ticks=false)
    hideydecorations!(axes.pressure_x; ticks=false)
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

    Colorbar(fig[2:4, 4], hts.advection_x; label=L"\text{Term} / f^2L_D")

    prettyrecord(n, fig, filename, frames; record_kw...)

    return fig
end