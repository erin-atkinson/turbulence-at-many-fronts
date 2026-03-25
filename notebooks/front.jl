@inline function interpolate(x, xs, fs)
    i1 = findlast(xs .<= x)
    i1 == nothing && return fs[end]
    i2 = i1 + 1
    i2 > length(fs) && return fs[i1]
    df = (x - xs[i1]) * (fs[i2] - fs[i1]) / (xs[i2] - xs[i1])
    return fs[i1] + df
end



function avg_buoyancy(DFM, z1, z2)
    k1, k2 = zᶜbounds(DFM, z1, z2)
    iterations, times = iterations_times(DFM)

    surface_mean(field) = mean(field[:, k1:k2]; dims=2)[:, 1]
    b_surface = timeseries_of(surface_mean, DFM, "b_coarse", iterations)

    return b_surface
end

function mixed_layer_height(DFM, b_surface, db)
    iterations, times = iterations_times(DFM)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(DFM)
    
    η_k = zeros(size(b_surface))
    η = zeros(size(b_surface))
    
    for n in axes(b_surface, 1)
        b = get_field(DFM, "b_dfm", iterations[n])
        for i in axes(b_surface, 2)
            b_depth = b[i, :] .- b_surface[n, i]
            η_k[n, i] = interpolate(db, b_depth, 1:length(zsᶜ))
            η[n, i] = interpolate(db, b_depth, zsᶜ)
        end
    end
    return η, η_k
end

function mld_buoyancy(DFM, η_k, z2)
    k2 = zᶜbounds(DFM, z2)
    iterations, times = iterations_times(DFM)
    
    b_mixed = zeros(size(η_k))
    
    for n in axes(η_k, 1)
        b = get_field(DFM, "b_dfm", iterations[n])
        for i in axes(η_k, 2)
            k0 = Int(floor(η_k[n, i]))
            k1 = k0 + 1
            b_mixed[n, i] = (sum(b[i, k1:k2]) + b[i, k0] * (η_k[n, i] - k0)) / (k2 - k0)
        end
    end

    return b_mixed
end

function cumulative_avg_b(file)
    iterations, times = iterations_times(file)
    sp = simulation_parameters(file)
    Δt = times[2] - times[1]
    
    return -cumsum([surface_b_flux(t, sp) * Δt / sp.H for t in times])
end

function front_width(file, b_front, bl, br)
    iterations, times = iterations_times(file)
    sp = simulation_parameters(file)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(file)
    
    ilr = zeros(size(b_front, 1), 2)
    xlr = zeros(size(b_front, 1), 2)
    ℓ = zeros(size(b_front, 1))
    
    for i in 1:length(times)
        il = interpolate(bl, b_front[i, :], 1:length(xsᶜ))
        ir = interpolate(br, b_front[i, :], 1:length(xsᶜ))
        
        xl = interpolate(bl, b_front[i, :], xsᶜ)
        xr = interpolate(br, b_front[i, :], xsᶜ)
        
        ilr[i, :] .= [il, ir]
        xlr[i, :] .= [xl, xr]
        ℓ[i] = xr - xl
    end
    return ilr, xlr, ℓ
end