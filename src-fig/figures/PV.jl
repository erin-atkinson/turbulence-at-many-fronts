# PV.jl


function interpolate_fluxes!(u, v, xs, ys, Jx, Jy)

    for i in axes(xs, 1), j in axes(ys, 2)
        u[i, j] = interpolate((xs[i], ys[j]), Jx)
        v[i, j] = interpolate((xs[i], ys[j]), Jy)
    end

    return nothing
end

@doc raw"""
    potential_vorticity_figure(run_id, frames;
        record_kw = NamedTuple(),
        N_window = 1,
        filename = joinpath(run_id, "$run_id-pv")
    )
Create a figure of the potential vorticity and its transport
"""
function potential_vorticity_figure(run_id, frames;
    record_kw = NamedTuple(),
    N_window = 1,
    filename = joinpath(run_id, "$run_id-pv")
    )
    
    PV = filepath(run_id, "PV", N_window)

    fts = fts_tuple(PV; q="q", Jx="Jx", Jz="Jz", 𝔍x="𝔍x", 𝔍z="𝔍z")

    sp = simulation_parameters(PV)
    times = fts.q.times
    
    n = Observable(frames[1])
    t = @lift interp_time($n, times)
    
    title = @lift let t_hr = @sprintf "%.0f" ($t / 3600)
        L"\text{Potential vorticity} \quad t = %$t_hr \, \text{hr}"
    end

    field_observables = make_fts_observables(fts, :, j, :, t; unit = sp.f * sp.N₀²)

    # Prepare fluxes
    xs_flux = range(-sp.Lh/2, sp.Lh/2, 32)
    zs_flux = range(-sp.Lz, 0, 16)
    u_advective = Observable([0.0 for x in xs_flux, z in zs_flux])
    w_advective = Observable([0.0 for x in xs_flux, z in zs_flux])
    u_diffusive = Observable([0.0 for x in xs_flux, z in zs_flux])
    w_diffusive = Observable([0.0 for x in xs_flux, z in zs_flux])

    on(n) do
        interpolate_fluxes!(u_advective, w_advective, xs_flux, zs_flux, Jx[], Jz[])
        interpolate_fluxes!(u_diffusive, w_diffusive, xs_flux, zs_flux, 𝔍x[], 𝔍z[])

        notify(u_advective)
        notify(w_advective)
        notify(u_diffusive)
        notify(w_diffusive)
    end


    fig = Figure(; size=(figure_width, 400), fontsize)
    Label(fig[1, 1:3], title)
    
    ax_kw = (;
        xlabel = x_label,
        ylabel = z_label,
        limits = transect_limits(sp)
    )

    ax_advective = Axis(fig[2, 1]; ax_kw...)
    ax_diffusive = Axis(fig[2, 2]; ax_kw...)

    hideydecorations!(ax_diffusive; ticks=false)
    
    ht_q = begin
        xs = nov(xnodes(fts.q; with_halos=true)) ./ x_unit
        zs = nov(znodes(fts.q; with_halos=true)) ./ z_unit
        data = field_observables.q
        colormap = :curl
        colorrange = (-0.25, 0.25)

        heatmap!(ax_advective, xs, zs, data; colormap, colorrange)
        heatmap!(ax_diffusive, xs, zs, data; colormap, colorrange)
    end

    begin
        xs = xs_flux ./ x_unit
        zs = zs_flux ./ z_unit
        
        arrows2d!(ax_advective, xs_flux, zs_flux, u_advective, w_advective)
        arrows2d!(ax_diffusive, xs_flux, zs_flux, u_diffusive, w_diffusive)
    end

    Colorbar(fig[2, 3], ht_q; label=q_label)

    colgap!(fig.layout, 40)
    prettyrecord(n, fig, filename, frames; record_kw...)

    return fig
end
@doc raw"""
    potential_vorticity_vertical_transport(run_id, frames, filename=nothing;
        record_kw = NamedTuple(),
        N_window = 1
    )
Create a 2D timeseries of potential vorticity vertical transport
"""
