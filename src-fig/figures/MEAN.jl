function mean_fields(foldername, frames, filename;
        fig_kw = NamedTuple(),
        ax_kw = NamedTuple(),
        record_kw = NamedTuple(),
        background = true)

        MEAN = joinpath(foldername, "MEAN.jld2")

        fds = FieldDataset(MEAN; backend=InMemory(10))
        sp = fds.metadata["parameters"]
        
        n = Observable(frames[1])

        title = 
        u_bar = @lift fds.u_bar[$n][:, 1, :] .* 100
        v_bar = @lift fds.v_bar[$n][:, 1, :] .* 100
        w_bar = @lift fds.w_bar[$n][:, 1, :] .* 1000
        b_bar = @lift fds.b_bar[$n][:, 1, :] ./ sp.Δb
        
        fig = Figure(; 
            size=(1000, 400),
            fig_kw...
        )

        ax_kw = (;
            xlabel = L"x / \text{km}",
            ylabel = L"z / \text{m}",
            limits = (-sp.Lh / 2000, sp.Lh / 2000, sp.Lz, 0)
        )

        ax_u = Axis(fig[2, 1]; aw_kw...)
        ax_v = Axis(fig[2, 2]; aw_kw...)
        ax_w = Axis(fig[2, 3]; aw_kw...)

        hideydecorations!(ax_v; ticks=false)
        hideydecorations!(ax_w; ticks=false)

        ht_u = begin
            xs = xnodes(fds.u_bar; with_halos=true)
            zs = znodes(fds.u_bar; with_halos=true)
            data = u_bar
            colormap = :balance
            colorrange = (-10, 10)

            heatmap!(ax_u, xs, zs, data; colormap, colorrange)
        end

        ht_v = begin
            xs = xnodes(fds.v_bar; with_halos=true)
            zs = znodes(fds.v_bar; with_halos=true)
            data = v_bar
            colormap = :balance
            colorrange = (-10, 10)

            heatmap!(ax_v, xs, zs, data; colormap, colorrange)
        end

        ht_w = begin
            xs = xnodes(fds.w_bar; with_halos=true)
            zs = znodes(fds.w_bar; with_halos=true)
            data = w_bar
            colormap = :balance
            colorrange = (-10, 10)

            heatmap!(ax_w, xs, zs, data; colormap, colorrange)
        end

        begin 
            xs = xnodes(fds.b_bar; with_halos=true)
            zs = znodes(fds.b_bar; with_halos=true)
            data = b_bar
            levels = b_levels
            color = (:black, 0.5)

            contour!(ax_u, xs, zs, data; levels, color)
            contour!(ax_v, xs, zs, data; levels, color)
            contour!(ax_w, xs, zs, data; levels, color)
        end

        Colorbar(fig[3, 1], ht_u; flipaxis=false, horizontal=true, label=u_bar_label)
        Colorbar(fig[3, 2], ht_v; flipaxis=false, horizontal=true, label=v_bar_label)
        Colorbar(fig[3, 3], ht_w; flipaxis=false, horizontal=true, label=w_bar_label)

end

function mean_hovmoller()

end