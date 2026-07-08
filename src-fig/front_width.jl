function front_width(xs, bs, bl, br)
    i_min = argmin(bs)
    b_front = bs .- minimum(bs)
    b_front = map(1:size(b_front, 1)) do i
        i > i_min ? b_front[i] : zero(eltype(b_front))
    end
    
    il = interpolate(bl, b_front, 1:length(xs))
    ir = interpolate(br, b_front, 1:length(xs))
    
    xl = interpolate(bl, b_front, xs)
    xr = interpolate(br, b_front, xs)
    
    ilr = [il, ir]
    xlr = [xl, xr]
    ℓ = xr - xl
    
    return ilr, xlr, ℓ
end

# Front width over time, averaged over depth
function average_front_width(COARSE)
    iterations, times = iterations_times(COARSE)
    sp = simulation_parameters(COARSE)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(COARSE)
    k1, k2 = zᶜbounds(COARSE, -sp.Lz, -5)
    bl = 1/12 * sp.Δb
    br = bl + 1/6 * sp.Δb
    
    N = length(iterations)
    
    ilr = zeros(N, 2)
    xlr = zeros(N, 2)
    ℓ = zeros(N)
    
    for i in 1:N
        b = get_field(COARSE, "b_coarse", iterations[i]) do field
            mean(field[:, k1:k2]; dims=2)[:, 1]
        end
        _ilr, _xlr, _ℓ = front_width(xsᶜ, b, bl, br)
        ilr[i, :] .= _ilr
        xlr[i, :] .= _xlr
        ℓ[i, :] .= _ℓ
    end

    return ilr, xlr, 6ℓ
end

function surface_front_width(COARSE)
    iterations, times = iterations_times(COARSE)
    sp = simulation_parameters(COARSE)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(COARSE)
    k1, k2 = zᶜbounds(COARSE, -10, -5)
    bl = 1/12 * sp.Δb
    br = bl + 1/6 * sp.Δb

    N = length(iterations)
    
    ilr = zeros(N, 2)
    xlr = zeros(N, 2)
    ℓ = zeros(N)
    
    for i in 1:N
        b = get_field(COARSE, "b_coarse", iterations[i]) do field
            mean(field[:, k1:k2]; dims=2)[:, 1]
        end
        _ilr, _xlr, _ℓ = front_width(xsᶜ, b, bl, br)
        ilr[i, :] .= _ilr
        xlr[i, :] .= _xlr
        ℓ[i, :] .= _ℓ
    end

    return ilr, xlr, 6ℓ
end

function mixed_front_width(MLD)
iterations, times = iterations_times(MLD)
    sp = simulation_parameters(MLD)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(MLD)
    bl = 1/12 * sp.Δb
    br = bl + 1/6 * sp.Δb
    
    N = length(iterations)
    
    ilr = zeros(N, 2)
    xlr = zeros(N, 2)
    ℓ = zeros(N)
    
    for i in 1:N
        b = get_field(MLD, "b_coarse_constant_mld", iterations[i])
        _ilr, _xlr, _ℓ = front_width(xsᶜ, b, bl, br)
        ilr[i, :] .= _ilr
        xlr[i, :] .= _xlr
        ℓ[i, :] .= _ℓ
    end

    return ilr, xlr, 6ℓ
end
