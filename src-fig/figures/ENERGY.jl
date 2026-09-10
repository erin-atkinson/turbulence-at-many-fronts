# ENERGY.jl

@doc raw"""
    energy_balance_check(run_id; N_window=1)
Return a figure that verifies the energy balance.

This function returns a figure that contains a timeseries for each term in the mean kinetic and potential energies
```math
\frac{\text{d}}{\text{d}t}\text{MKE} = \text{DSP} + \text{WIND} + \text{BUOYANCY} + \text{SPONGE}_\text{MKE} - \text{LSP} - \text{VSP} + \text{STRAIN}_\text{MKE}
```

```math
\frac{\text{d}}{\text{d}t}\text{MPE} = -\text{BUOYANCY} + \text{SPONGE}_\text{MPE} - \text{BFLUX} + \text{MIXED} + \text{COOLING} + \text{STRAIN}_\text{MPE}
```
"""
function energy_balance_check(run_id; N_window=1)
    mkebalance = MKEBalance(run_id, N_window)
    mpebalance = MPEBalance(run_id, N_window)

    _, times = iterations_times(mkebalance)

    tendency_unit = 1/1037

    fig = Figure(; size=(figure_width, 300), fontsize)

    ax = Axis(fig[1, 1]; 
        xlabel = t_label,
        ylabel = L"\text{kW} \, \text{km}^{-1}",
        limits = (0, times[end] / t_unit, nothing, nothing),
    )

    lns = plot_balance!(ax, times ./ t_unit, mkebalance, tendency_unit)
    plot_target!(ax, times ./ t_unit, mkebalance, tendency_unit)
    plot_total!(ax, times ./ t_unit, mkebalance, tendency_unit)
    make_legend!(fig[1, 2], lns, mkebalance; title=L"A") 
    
    ax = Axis(fig[2, 1]; 
        xlabel = t_label,
        ylabel = L"\text{kW} \, \text{km}^{-1}",
        limits = (0, times[end] / t_unit, nothing, nothing),
    )

    lns = plot_balance!(ax, times ./ t_unit, mpebalance, tendency_unit)
    plot_target!(ax, times ./ t_unit, mpebalance, tendency_unit)
    plot_total!(ax, times ./ t_unit, mpebalance, tendency_unit)
    make_legend!(fig[2, 2], lns, mpebalance; title=L"A") 

    fig
end

function terms_MKE(run_id, frames, filename;
    fig_kw = NamedTuple(),
    ax_kw = NamedTuple(),
    record_kw = NamedTuple(),
    N_window = 1
    )
    foldername = joinpath(scratchpath, run_id)
    
    suffix = N_window == 1 ? "" : "-$N_window"
    MEAN = joinpath(foldername, "MEAN$(suffix).jld2")
    ENERGY = joinpath(foldername, "ENERGY$(suffix).jld2")
    
    sp = simulation_parameters(MEAN)
    iterations, times = iterations_times(MEAN)
    
    n = Observable(frames[1])
    t = @lift interp_time($n, times)
    
    fts_b_bar = FieldTimeSeries(MEAN, "b_bar")
    
    fieldnames = MKE_density_terms

    fts = NamedTuple(k => FieldTimeSeries(ENERGY, v; backend=OnDisk()) for (k, v) in pairs(fieldnames))
    
    title = @lift let t_hr = @sprintf "%.0f" ($t / 3600)
        L"\text{Terms in MKE balance}\quad t = %$t_hr \, \text{hr}"
    end
    
    b_bar = @lift nov(fts_b_bar[Time($t)][:, 1, :]) ./ sp.Δb

    field_data = NamedTuple(k => @lift MKE_signs[k] * nov(v[Time($t)][:, 1, :]) ./ (sp.L^2 * sp.f^3) for (k, v) in pairs(fts))
    
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
        dsp = Axis(fig[2, 1]; ax_kw..., title=MKE_term_labels.dsp),
        lsp = Axis(fig[2, 2]; ax_kw..., title=MKE_term_labels.lsp),
        vsp = Axis(fig[2, 3]; ax_kw..., title=MKE_term_labels.vsp),
        buoyancy = Axis(fig[3, 1]; ax_kw..., title=MKE_term_labels.buoyancy),
        sponge_mke = Axis(fig[3, 2]; ax_kw..., title=MKE_term_labels.sponge_mke),
        strain_mke = Axis(fig[3, 3]; ax_kw..., title=MKE_term_labels.strain_mke),
    )

    hidexdecorations!(axes.dsp; ticks=false)
    hidexdecorations!(axes.lsp; ticks=false)
    hidexdecorations!(axes.vsp; ticks=false)

    hideydecorations!(axes.lsp; ticks=false)
    hideydecorations!(axes.vsp; ticks=false)
    hideydecorations!(axes.sponge_mke; ticks=false)
    hideydecorations!(axes.strain_mke; ticks=false)

    hts = NamedTuple(
        map(keys(axes)) do k
            xs = nov(xnodes(fts[k]; with_halos=true)) ./ x_unit
            zs = nov(znodes(fts[k]; with_halos=true)) ./ z_unit
            data = field_data[k]
            colormap = :balance
            colorrange = (-0.1, 0.1)
    
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

    Colorbar(fig[2:3, 4], hts.dsp; label=L"\text{Term} / f^3L_D^2")

    prettyrecord(n, fig, filename, frames; record_kw...)

    return fig
end