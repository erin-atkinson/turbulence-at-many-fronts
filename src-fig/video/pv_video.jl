function pv_video(foldername, frames, filename;
        fig_kw = NamedTuple(),
        ax_kw = NamedTuple(),
        ht_kw = NamedTuple(),
        record_kw = NamedTuple(),
        σ = 0
    )

    PV = jldopen(joinpath(foldername, "PV.jld2"))

    iterations, times = iterations_times(PV)
    ∫αdt = normalized_time(PV)
    sp = simulation_parameters(PV)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(PV)
    inds = center_indices(PV)

    n = Observable(frames[1])
    t = @lift times[$n]
    iteration = @lift iterations[$n]
    
    fig = Figure(; size=(800, 540), fontsize=18)

    C = @sprintf "%02.1f" :C in keys(sp) ? sp.C : sp.Ek
    Roa = @sprintf "%02.2f" sp.Roα

    title = @lift let
        hr_val = @sprintf "%03.0f" times[$n] / 3600
        s_val = @sprintf "%03.2f" ∫αdt[$n]
        L"(%$Roa, %$C) \quad s = %$(s_val) \quad t = %$(hr_val)\,\text{hr}"
    end

    q = @lift get_field(PV, "coarse_q", iterations[$n]) do field
        filt(field, σ) ./ (sp.N₀² * sp.f)
    end

    ax = Axis(fig[1, 1];
        xlabel = L"x / \text{km}",
        ylabel = L"z / \text{m}",
        limits = (-sp.Lh / 2000, sp.Lh / 2000, -sp.Lz, 0),
        title,
        xticks = [-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5]
    )

    colormap = to_colormap(:curl)
    
    ht = heatmap!(ax, xsᶠ ./ 1000, zsᶠ, q; 
        colormap,
        highclip = colormap[end],
        lowclip = colormap[1],
        colorrange = (-0.1, 0.1),
        ht_kw...
    )
    
    Colorbar(fig[1, 2], ht; label=L"q / N_0^2f")
    
    length(frames) > 1 && prettyrecord(n, fig, filename, frames; record_kw...)
    return fig
end