function uv_video(foldername, frames, filename, z;
        fig_kw = NamedTuple(),
        ax_kw = NamedTuple(),
        axh_kw = NamedTuple(),
        htu_kw = NamedTuple(),
        htv_kw = NamedTuple(),
        ct_kw = NamedTuple(),
        record_kw = NamedTuple(),
        background = true,
        coarse = false
    )

    file = joinpath(foldername, coarse ? "COARSE.jld2" : "DFM.jld2")
    FILE = jldopen(file)
    OUTPUT = jldopen(joinpath(foldername, "output.jld2"))

    iterations, times = iterations_times(FILE)
    ∫αdt = normalized_time(FILE)
    sp = simulation_parameters(FILE)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(FILE)
    inds = center_indices(FILE)
    k = zᶜbounds(FILE, z)

    n = Observable(frames[1])
    t = @lift times[$n]
    iteration = @lift iterations[$n]
    
    fig = Figure(; size=(1000, 540), fontsize=18)

    C = @sprintf "%02.1f" :C in keys(sp) ? sp.C : sp.Ek
    Roa = @sprintf "%02.2f" sp.Roα

    title = @lift let
        hr_val = @sprintf "%03.0f" times[$n] / 3600
        s_val = @sprintf "%03.2f" ∫αdt[$n]
        L"(%$Roa, %$C) \quad s = %$(s_val) \quad t = %$(hr_val)\,\text{hr}"
    end

    u_name = coarse ? "u_coarse" : "u_dfm"
    v_name = coarse ? "v_coarse" : "v_dfm"
    b_name = coarse ? "b_coarse" : "b_dfm"

    U = @lift background .* [-100variable_strain_rate($t, sp) * x for x in xsᶠ, y in 1:1]

    u = @lift get_field(FILE, u_name, $iteration) do field
        field .* 100 .+ $U
    end
    uh = @lift get_field(OUTPUT, "u", $iteration) do field
        field[:, :, k] .* 100 .+ $U
    end

    v = @lift get_field(FILE, v_name, $iteration) do field
        field .* 100
    end
    vh = @lift get_field(OUTPUT, "v", $iteration) do field
        field[:, :, k] .* 100
    end
    
    b = @lift get_field(FILE, b_name, $iteration)
    bh = @lift get_field(OUTPUT, "b", $iteration) do field
        field[:, :, k]
    end

    levels = @lift begin
        bmin = min(minimum($b), minimum($bh))
        bmax = max(maximum($b), maximum($bh))
        bmin:(sp.Δb / 6):bmax
    end
    
    ax_kw = (; 
        xlabel = L"x / \text{km}",
        ylabel = L"z / \text{m}",
        limits = (-sp.Lh / 2000, sp.Lh / 2000, -sp.Lz, 0),
        title,
        xticks = [-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5],
        ax_kw...
    )

    axh_kw = (; 
        xlabel = L"x / \text{km}",
        ylabel = L"y / \text{km}",
        limits = (-sp.Lh / 2000, sp.Lh / 2000, -sp.Ly/2000, sp.Ly/2000),
        xticks = [-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5],
        yticks = [-0.2, -0.1, 0, 0.1, 0.2],
        axh_kw...
    )

    ax_u = Axis(fig[1, 1]; ax_kw...)
    ax_v = Axis(fig[2, 1]; ax_kw..., title="")
    hidexdecorations!(ax_u)

    axh_u = Axis(fig[1, 2]; axh_kw...)
    axh_v = Axis(fig[2, 2]; axh_kw..., title="")
    hidexdecorations!(axh_u)
    
    ct_kw = (; 
        levels,
        color = (:black, 0.5),
        ct_kw...
    )
    colormap = to_colormap(:balance)
    htu_kw = (; 
        colormap,
        colorrange = (-5, 5),
        highclip = colormap[end],
        lowclip = colormap[1],
    )

    htv_kw = (; 
        colormap,
        colorrange = (-10, 10),
        highclip = colormap[end],
        lowclip = colormap[1],
    )

    heatmap!(ax_u, xsᶠ ./ 1000, zsᶜ, u; htu_kw...)
    heatmap!(ax_v, xsᶜ ./ 1000, zsᶜ, v; htv_kw...)

    htu = heatmap!(axh_u, xsᶠ ./ 1000, ysᶜ ./ 1000, uh; htu_kw...)
    htv = heatmap!(axh_v, xsᶜ ./ 1000, ysᶠ ./ 1000, vh; htv_kw...)
    
    contour!(ax_u, xsᶜ ./ 1000, zsᶜ, b; ct_kw...)
    contour!(ax_v, xsᶜ ./ 1000, zsᶜ, b; ct_kw...)

    #contour!(axh_u, xsᶜ ./ 1000, ysᶜ ./ 1000, bh; ct_kw...)
    #contour!(axh_v, xsᶜ ./ 1000, ysᶜ ./ 1000, bh; ct_kw...)
    
    Colorbar(fig[1, 3], htu; label=L"u / \text{cm}\,\text{s}^{-1}")
    Colorbar(fig[2, 3], htv; label=L"v / \text{cm}\,\text{s}^{-1}")

    subfig_label!(fig[1, 1], 1; fontsize=18)
    subfig_label!(fig[1, 2], 2; fontsize=18)
    subfig_label!(fig[2, 1], 3; fontsize=18)
    subfig_label!(fig[2, 2], 4; fontsize=18)
    
    length(frames) > 1 && prettyrecord(n, fig, filename, frames; record_kw...)
    return fig
end