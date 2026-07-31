MKE_term_labels = (;
    dsp = L"\text{DSP}",
    lsp = L"-\text{LSP}",
    vsp = L"-\text{VSP}",
    buoyancy = L"\text{BUOYANCY}",
    sponge_mke = L"\text{SPONGE}_\text{MKE}",
    strain_mke = L"\text{STRAIN}_\text{MKE}",
    wind = L"\text{WIND}"
)

MKE_terms = (;
    dsp = "dsp",
    lsp = "lsp",
    vsp = "vsp",
    buoyancy = "buoyancy",
    sponge_mke = "sponge_mke",
    strain_mke = "strain_mke",
    wind = "wind"
)

MKE_density_terms = (;
    dsp = "dsp_density",
    lsp = "lsp_density",
    vsp = "vsp_density",
    buoyancy = "buoyancy_density",
    sponge_mke = "sponge_mke_density",
    strain_mke = "strain_mke_density"
)

MKE_signs = (;
    dsp = 1,
    lsp = -1,
    vsp = -1,
    buoyancy = 1,
    sponge_mke = 1,
    strain_mke = 1,
    wind = 1,
)

MPE_term_labels = (;
    buoyancy = L"-\text{BUOYANCY}",
    bflux = L"-\text{BFLUX}",
    cooling = L"\text{COOLING}",
    strain_mpe = L"\text{STRAIN}_\text{MPE}",
    mixed = L"\text{MIXED}",
)

MPE_terms = (;
    buoyancy = "buoyancy",
    bflux = "bflux",
    cooling = "cooling",
    strain_mpe = "strain_mpe",
    mixed = "mixed"
)

MPE_signs = (;
    buoyancy = -1,
    bflux = -1,
    cooling = 1,
    strain_mpe = 1,
    mixed = 1
)

function check_ENERGY(run_id; N_window=1)
    foldername = joinpath(scratchpath, run_id)
    
    suffix = N_window == 1 ? "" : "-$N_window"
    MEAN = joinpath(foldername, "MEAN$(suffix).jld2")
    ENERGY = joinpath(foldername, "ENERGY$(suffix).jld2")
        
    iterations, times = iterations_times(ENERGY)
    sp = simulation_parameters(ENERGY)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(ENERGY)
    Δt = times[2] - times[1]
    fig = Figure(; size=(600, 600), fontsize=18)
    
    ax_mke = Axis(fig[1, 1]; 
        xlabel = t_label,
        ylabel = L"\text{kW} \, \text{km}^{-1}",
        limits = (0, times[end] / t_unit, nothing, nothing)
    )

    ax_mpe = Axis(fig[2, 1]; 
        xlabel = t_label,
        ylabel = L"\text{kW} \, \text{km}^{-1}",
        limits = (0, times[end] / t_unit, nothing, nothing)
    )

    energy_unit = 1/1037 #sp.L^2 * sp.f^2 * sp.L * sp.H
    power_unit = energy_unit / 1
    
    mke = timeseries_of(ENERGY, "mke", iterations)
    mpe = timeseries_of(ENERGY, "mpe", iterations)
    
    total_mke_actual = diff(mke) ./ Δt ./ power_unit
    total_mpe_actual = diff(mpe) ./ Δt ./ power_unit

    mke_terms = NamedTuple(k => MKE_signs[k] * timeseries_of(ENERGY, MKE_terms[k], iterations) ./ power_unit for k in keys(MKE_terms))
    mpe_terms = NamedTuple(k => MPE_signs[k] * timeseries_of(ENERGY, MPE_terms[k], iterations) ./ power_unit for k in keys(MPE_terms))

    mke_lines = NamedTuple(k => lines!(ax_mke, times ./ t_unit, mke_terms[k]) for k in keys(mke_terms))
    mpe_lines = NamedTuple(k => lines!(ax_mpe, times ./ t_unit, mpe_terms[k]) for k in keys(mpe_terms))

    total_mke = sum(mke_terms)
    total_mpe = sum(mpe_terms)
    
    lines!(ax_mke, times ./ t_unit, total_mke; color=:black)
    lines!(ax_mke, times[2:end] ./ t_unit, total_mke_actual; color=:black, linestyle=:dash)
    lines!(ax_mpe, times ./ t_unit, total_mpe; color=:black)
    lines!(ax_mpe, times[2:end] ./ t_unit, total_mpe_actual; color=:black, linestyle=:dash)

    Legend(fig[1, 2], [ln for ln in mke_lines], [ln for ln in MKE_term_labels])
    Legend(fig[2, 2], [ln for ln in mpe_lines], [ln for ln in MPE_term_labels])
    
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