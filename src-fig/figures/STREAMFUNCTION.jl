STREAMFUNCTION_term_labels = (;
    ageostrophic = L"\text{Ageostrophic}",
    strain = L"\text{Strain}",
    turbulence = L"\text{Turbulence}",
    sponge = L"\text{Sponge}",
    sce = L"\text{SCE}",
)

# Terms in the balance equation only 
STREAMFUNCTION_terms = (;
    ageostrophic = "ageostrophic",
    strain = "strain",
    turbulence = "turbulence",
    sponge = "sponge",
)

function check_STREAMFUNCTION(run_id; N_window=1)
    foldername = joinpath(scratchpath, run_id)
    
    suffix = N_window == 1 ? "" : "-$N_window"
    MEAN = joinpath(foldername, "MEAN$(suffix).jld2")
    STREAMFUNCTION = joinpath(foldername, "STREAMFUNCTION$(suffix).jld2")
        
    iterations, times = iterations_times(STREAMFUNCTION)
    sp = simulation_parameters(STREAMFUNCTION)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(STREAMFUNCTION)
    Δt = times[2] - times[1]
    fig = Figure(; size=(600, 400), fontsize=18)
    
    ax_sce = Axis(fig[1, 1]; 
        xlabel = t_label,
        ylabel = L"\text{kW} \, \text{km}^{-1}",
        limits = (0, times[end] / t_unit, nothing, nothing)
    )

    energy_unit = 1/1037 #sp.L^2 * sp.f^2 * sp.L * sp.H
    power_unit = energy_unit / 1
    
    sce = timeseries_of(STREAMFUNCTION, "sce", iterations)
    
    total_sce_actual = diff(sce) ./ Δt ./ power_unit

    sce_terms = NamedTuple(k => timeseries_of(STREAMFUNCTION, STREAMFUNCTION_terms[k], iterations) ./ power_unit for k in keys(STREAMFUNCTION_terms))

    sce_lines = NamedTuple(k => lines!(ax_sce, times ./ t_unit, sce_terms[k]) for k in keys(sce_terms))

    total_sce = sum(sce_terms)
    
    lines!(ax_sce, times ./ t_unit, total_sce; color=:black)
    lines!(ax_sce, times[2:end] ./ t_unit, total_sce_actual; color=:black, linestyle=:dash)

    Legend(fig[1, 2], [ln for ln in sce_lines], [STREAMFUNCTION_term_labels[k] for k in keys(sce_terms)])
    
    fig
end

function compare_STREAMFUNCTION(run_ids; N_window=1, normalize=false)
    suffix = N_window == 1 ? "" : "-$N_window"
    fig = Figure(; fontsize=18, size=(800, 600))

    energy_label = normalize ? L"BHL / f" : L"\text{kJ} \, \text{km}^{-1}"
    power_label = normalize ? L"BHL" : L"\text{kW} \, \text{km}^{-1}"
    
    axes = (;
        sce = Axis(fig[1, 1]; xlabel=t_label, ylabel=energy_label, title=STREAMFUNCTION_term_labels.sce),
        ageostrophic = Axis(fig[1, 2]; xlabel=t_label, ylabel=power_label, title=STREAMFUNCTION_term_labels.ageostrophic),
        strain = Axis(fig[1, 3]; xlabel=t_label, ylabel=power_label, title=STREAMFUNCTION_term_labels.strain),
        turbulence = Axis(fig[2, 1]; xlabel=t_label, ylabel=power_label, title=STREAMFUNCTION_term_labels.turbulence),
        sponge = Axis(fig[2, 2]; xlabel=t_label, ylabel=power_label, title=STREAMFUNCTION_term_labels.sponge),
    )
    
    hidexdecorations!(axes.sce; ticks=false, grid=false)
    hidexdecorations!(axes.ageostrophic; ticks=false, grid=false)
    
    lns = map(run_ids) do run_id
        foldername = joinpath(scratchpath, run_id)
        STREAMFUNCTION = joinpath(foldername, "STREAMFUNCTION$suffix.jld2")
        
        iterations, times = iterations_times(STREAMFUNCTION)
        sp = simulation_parameters(STREAMFUNCTION)

        energy_unit = normalize ? sp.B * sp.H * sp.L / sp.f : 1/1037
        power_unit = normalize ? sp.B * sp.H * sp.L : energy_unit / 1
        
        timeseries = timeseries_of(STREAMFUNCTION, iterations; STREAMFUNCTION_terms...)
        sce = timeseries_of(STREAMFUNCTION, "sce", iterations)
        
        map(keys(timeseries)) do k
            lines!(axes[k], times ./ t_unit, timeseries[k] ./ power_unit)
        end
        lines!(axes.sce, times ./ t_unit, sce ./ energy_unit)
    end

    Legend(fig[2, 3], lns, run_ids; tellwidth=false)
    fig
end

function terms_STREAMFUNCTION(run_id, frames, filename;
    fig_kw = NamedTuple(),
    ax_kw = NamedTuple(),
    ax_z_kw = NamedTuple(),
    record_kw = NamedTuple(),
    N_window = 1
    )
    foldername = joinpath(scratchpath, run_id)
    
    suffix = N_window == 1 ? "" : "-$N_window"
    STREAMFUNCTION = joinpath(foldername, "STREAMFUNCTION$suffix.jld2")
    MEAN = joinpath(foldername, "MEAN$suffix.jld2")

    sp = simulation_parameters(MEAN)
    iterations, times = iterations_times(MEAN)
    
    n = Observable(frames[1])
    t = @lift interp_time($n, times)
    
    fts_b_bar = FieldTimeSeries(MEAN, "b_bar")
    
    fieldnames = (;
        advection_self = "advection_self",
        advection_background = "advection_background",
        ageostrophic = "ageostrophic",
        strain = "strain",
        turbulence = "turbulence",
        sponge = "sponge"
    )

    fts = NamedTuple(k => FieldTimeSeries(STREAMFUNCTION, v; backend=OnDisk()) for (k, v) in pairs(fieldnames))

    title = @lift let t_hr = @sprintf "%.0f" ($t / 3600)
        L"\text{Terms in }\nabla^2 \psi\text{ balance}\quad t = %$t_hr \, \text{hr}"
    end
    
    b_bar = @lift nov(fts_b_bar[Time($t)][:, 1, :]) ./ sp.Δb

    field_data = NamedTuple(k => @lift nov(v[Time($t)][:, 1, :]) ./ sp.f^2 for (k, v) in pairs(fts))
    
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
        advection_self = Axis(fig[2, 1]; ax_kw..., title=STREAMFUNCTION_term_labels.advection_self),
        advection_background = Axis(fig[2, 2]; ax_kw..., title=STREAMFUNCTION_term_labels.advection_background),
        ageostrophic = Axis(fig[2, 3]; ax_kw..., title=STREAMFUNCTION_term_labels.ageostrophic),
        strain = Axis(fig[3, 1]; ax_kw..., title=STREAMFUNCTION_term_labels.strain),
        turbulence = Axis(fig[3, 2]; ax_kw..., title=STREAMFUNCTION_term_labels.turbulence),
        sponge = Axis(fig[3, 3]; ax_kw..., title=STREAMFUNCTION_term_labels.sponge)
    )

    hidexdecorations!(axes.advection_self; ticks=false)
    hidexdecorations!(axes.advection_background; ticks=false)
    hidexdecorations!(axes.ageostrophic; ticks=false)

    hideydecorations!(axes.advection_background; ticks=false)
    hideydecorations!(axes.ageostrophic; ticks=false)
    hideydecorations!(axes.turbulence; ticks=false)
    hideydecorations!(axes.sponge; ticks=false)

    
    hts = NamedTuple(
        map(keys(axes)) do k
            xs = nov(xnodes(fts[k]; with_halos=true)) ./ x_unit
            zs = nov(znodes(fts[k]; with_halos=true)) ./ z_unit
            data = field_data[k]
            colormap = :balance
            colorrange = (-10, 10)
    
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

    Colorbar(fig[2:3, 4], hts.advection_self; label=L"\text{Term} / f^2")

    #colgap!(fig.layout, 40)
    prettyrecord(n, fig, filename, frames; record_kw...)

    return fig
end