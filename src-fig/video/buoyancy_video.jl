function buoyancy_video(foldername, frames, filename;
        fig_kw = NamedTuple(),
        ax_kw = NamedTuple(),
        ct_kw = NamedTuple(),
        record_kw = NamedTuple(),
        σ = 0
    )

    DFM = jldopen(joinpath(foldername, "DFM.jld2"))

    iterations, times = iterations_times(DFM)
    ∫αdt = normalized_time(DFM)
    sp = simulation_parameters(DFM)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(DFM)
    inds = center_indices(DFM)
    k1, k2 = zᶜbounds(DFM, -0.9sp.H, -0.1sp.H)

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

    b = @lift get_field(DFM, "b_dfm", iterations[$n]) do field
        (field .- mean(field[inds, k1:k2])) ./ sp.Δb
    end

    ax = Axis(fig[1, 1];
        xlabel = L"x / \text{km}",
        ylabel = L"z / \text{m}",
        limits = (-sp.Lh / 2000, sp.Lh / 2000, -sp.Lz, 0),
        title,
        xticks = [-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5]
    )

    levels = (-3:3) ./ 6
    colormap = reverse(to_colormap(:deep))
    
    ctf = contourf!(ax, xsᶜ ./ 1000, zsᶜ, b; 
        colormap,
        levels,
        extendlow=colormap[1],
        extendhigh=colormap[end],
        ct_kw...
    )
    
    Colorbar(fig[1, 2], ctf; label=L"b / \Delta b")
    
    length(frames) > 1 && prettyrecord(n, fig, filename, frames; record_kw...)
    return fig
end