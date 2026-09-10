# BBALANCE.jl

@doc raw"""
    b_balance_check(run_id; N_window=1)
Return a figure that verifies the buoyancy balance.

This function returns a figure that contains a timeseries for each term in the quadratic balance for the buoyancy
```math
\frac{\overline{\text{D}}_\alpha\overline b}{\text{D}t}\overline{b} = \left (\mathscr{D} -B\delta (z)\right ) \overline b
```
"""
function b_balance_check(run_id; N_window=1)
    balance = BBalance(run_id, N_window)

    _, times = iterations_times(balance)
    
    tendency_unit = 0.01^2
    fig = Figure(; size=(figure_width, 300), fontsize)
    ax = Axis(fig[1, 1]; 
        xlabel = t_label,
        ylabel = L"A / \text{cm}^{2} \, \text{s}^{-4}",
        limits = (0, times[end] / t_unit, nothing, nothing),
    )

    lns = plot_balance!(ax, times ./ t_unit, balance, tendency_unit)
    plot_target!(ax, times ./ t_unit, balance, tendency_unit)
    plot_total!(ax, times ./ t_unit, balance, tendency_unit)

    make_legend!(fig[1, 2], lns, balance; title=L"A") 

    fig
end

function terms_BBALANCE(run_id, frames, filename;
    fig_kw = NamedTuple(),
    ax_kw = NamedTuple(),
    ax_z_kw = NamedTuple(),
    record_kw = NamedTuple(),
    N_window = 1
    )
    foldername = joinpath(scratchpath, run_id)
    
    suffix = N_window == 1 ? "" : "-$N_window"
    MEAN = joinpath(foldername, "MEAN$suffix.jld2")
    BBALANCE = joinpath(foldername, "BBALANCE$suffix.jld2")

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
        mixing_z = "mixing_z"
    )

    fts = NamedTuple(k => FieldTimeSeries(BBALANCE, v; backend=OnDisk()) for (k, v) in pairs(fieldnames))

    # Also vertical profiles
    fts_flux_density_z = FieldTimeSeries(BBALANCE, "flux_density_z"; backend=OnDisk())
    fts_turbulent_flux_density_z = FieldTimeSeries(BBALANCE, "turbulent_flux_density_z"; backend=OnDisk())

    # Need to take an average
    temp_flux_density_z = similar(fts_flux_density_z[Time(t[])])
    temp_turbulent_flux_density_z = similar(fts_turbulent_flux_density_z[Time(t[])])
    
    condition(i, j, k, grid, c) = -sp.Lh < Oceananigans.Grids.xnode(i, grid, Center()) < sp.Lh
    mean_flux_density_z = Field(Average(temp_flux_density_z; dims=(1, 2), condition))
    mean_turbulent_flux_density_z = Field(Average(temp_turbulent_flux_density_z; dims=(1, 2), condition))

    flux_density_z = Observable(nov(mean_flux_density_z[1, 1, :]) ./ sp.B)
    turbulent_flux_density_z = Observable(nov(mean_turbulent_flux_density_z[1, 1, :]) ./ sp.B)
    
    on(t) do t
        set!(temp_flux_density_z, fts_flux_density_z[Time(t)])
        set!(temp_turbulent_flux_density_z, fts_turbulent_flux_density_z[Time(t)])
        
        compute!(mean_flux_density_z)
        compute!(mean_turbulent_flux_density_z)
        
        flux_density_z[] = nov(mean_flux_density_z[1, 1, :]) ./ sp.B
        turbulent_flux_density_z[] = nov(mean_turbulent_flux_density_z[1, 1, :]) ./ sp.B
    end
    
    title = @lift let t_hr = @sprintf "%.0f" ($t / 3600)
        L"\text{Terms in }\overline{b}\text{ balance}\quad t = %$t_hr \, \text{hr}"
    end
    
    b_bar = @lift nov(fts_b_bar[Time($t)][:, 1, :]) ./ sp.Δb

    field_data = NamedTuple(k => @lift nov(v[Time($t)][:, 1, :]) ./ (sp.Δb * sp.f) for (k, v) in pairs(fts))
    
    fig = Figure(; 
        size=(1000, 600),
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
        advection_x = Axis(fig[2, 1]; ax_kw..., title=BBALANCE_term_labels.advection_x),
        advection_background = Axis(fig[2, 2]; ax_kw..., title=BBALANCE_term_labels.advection_background),
        advection_z = Axis(fig[2, 3]; ax_kw..., title=BBALANCE_term_labels.advection_z),
        mixing_x = Axis(fig[3, 1]; ax_kw..., title=BBALANCE_term_labels.mixing_x),
        mixing_z = Axis(fig[3, 2]; ax_kw..., title=BBALANCE_term_labels.mixing_z)
    )
    
    ax_z_kw = (;
        xlabel = L"\text{Flux density in front} / B",
        ylabel = z_label,
        limits = (-1, 1, -sp.Lz / z_unit, 0),
    )
    
    ax_z = Axis(fig[3, 3]; ax_z_kw...)

    hidexdecorations!(axes.advection_x; ticks=false)
    hidexdecorations!(axes.advection_background; ticks=false)
    hidexdecorations!(axes.advection_z; ticks=false)

    hideydecorations!(axes.advection_background; ticks=false)
    hideydecorations!(axes.advection_z; ticks=false)
    hideydecorations!(axes.mixing_z; ticks=false)
    hideydecorations!(ax_z; ticks=false, grid=false)

    
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
    
    begin
        zs = nov(znodes(fts_flux_density_z; with_halos=true)) ./ z_unit
        data = flux_density_z
        lines!(ax_z, data, zs)
    end
    
    begin
        zs = nov(znodes(fts_turbulent_flux_density_z; with_halos=true)) ./ z_unit
        data = turbulent_flux_density_z
        lines!(ax_z, data, zs)
    end
    
    map(keys(axes)) do k
        xs = nov(xnodes(fts_b_bar; with_halos=true)) ./ x_unit
        zs = nov(znodes(fts_b_bar; with_halos=true)) ./ z_unit
        data = b_bar
        levels = b_levels(fts_b_bar, sp) ./ sp.Δb
        color = (:black, 0.5)

        contour!(axes[k], xs, zs, data; levels, color)
    end

    Colorbar(fig[2:3, 4], hts.advection_x; label=L"\text{Term} / f\Delta b")

    #colgap!(fig.layout, 40)
    prettyrecord(n, fig, filename, frames; record_kw...)

    return fig
end