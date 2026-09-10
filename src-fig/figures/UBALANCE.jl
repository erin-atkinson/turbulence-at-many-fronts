# UBALANCE.jl

@doc raw"""
    u_balance_check(run_id; N_window=1)
Return a figure that verifies the across-front velocity balance.

This function returns a figure that contains a timeseries for each term in the quadratic balance for the across-front velocity
```math
\frac{\overline{\text{D}}_\alpha\overline u}{\text{D}t}\overline{u} = \left (f\overline v-\frac{\partial \overline{p}}{\partial x} -
\frac{\partial U}{\partial x}\overline u+\mathscr{F}_u + \overline S_u + \tau_x\delta (z)\right ) \overline u
```
"""
function u_balance_check(run_id; N_window=1)
    balance = UBalance(run_id, N_window)

    _, times = iterations_times(balance)
    
    tendency_unit = 0.01^2
    fig = Figure(; size=(figure_width, 300), fontsize)
    ax = Axis(fig[1, 1]; 
        xlabel = t_label,
        ylabel = L"A / \text{cm}^{2} \, \text{s}^{-3}",
        limits = (0, times[end] / t_unit, nothing, nothing),
    )

    lns = plot_balance!(ax, times ./ t_unit, balance, tendency_unit)
    plot_target!(ax, times ./ t_unit, balance, tendency_unit)
    plot_total!(ax, times ./ t_unit, balance, tendency_unit)

    make_legend!(fig[1, 2], lns, balance; title=L"A") 

    fig
end

@doc raw"""
    u_balance_profiles(run_id, il, ir; N_window=1)
Return a figure that verifies the across-front velocity balance.

This function returns a figure that contains a timeseries for each term in the quadratic balance for the across-front velocity
```math
\frac{\overline{\text{D}}_\alpha\overline u}{\text{D}t}\overline{u} = \left (f\overline v-\frac{\partial \overline{p}}{\partial x} -
\frac{\partial U}{\partial x}\overline u+\mathscr{F}_u + \overline S_u + \tau_x\delta (z)\right ) \overline u
```
"""

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