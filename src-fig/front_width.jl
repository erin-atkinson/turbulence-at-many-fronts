# Frontal width averaged between z1, z2
function front_width(file, field, iterations, z1, z2, bl, br)
    iterations, times = iterations_times(file)
    sp = simulation_parameters(file)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(file)
    
    N = length(iterations)
    ilr = zeros(N, 2)
    xlr = zeros(N, 2)
    ℓ = zeros(N)

    k1, k2 = zᶜbounds(file, z1, z2)
    
    for i in 1:N
        b = get_field(file, field, iterations[i])
        b_front = mean(b[:, k1:k2]; dims=2)[:, 1]
        b_front .-= mean(b_front)
        
        il = interpolate(bl, b_front, 1:length(xsᶜ))
        ir = interpolate(br, b_front, 1:length(xsᶜ))
        
        xl = interpolate(bl, b_front, xsᶜ)
        xr = interpolate(br, b_front, xsᶜ)
        
        ilr[i, :] .= [il, ir]
        xlr[i, :] .= [xl, xr]
        ℓ[i] = xr - xl
    end
    
    return ilr, xlr, ℓ
end

function front_width(file, b, bl, br)
    iterations, times = iterations_times(file)
    sp = simulation_parameters(file)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(file)
    
    N = length(iterations)
    ilr = zeros(N, 2)
    xlr = zeros(N, 2)
    ℓ = zeros(N)
    
    for i in 1:N
        b_front = b[i, :]
        
        il = interpolate(bl, b_front, 1:length(xsᶜ))
        ir = interpolate(br, b_front, 1:length(xsᶜ))
        
        xl = interpolate(bl, b_front, xsᶜ)
        xr = interpolate(br, b_front, xsᶜ)
        
        ilr[i, :] .= [il, ir]
        xlr[i, :] .= [xl, xr]
        ℓ[i] = xr - xl
    end
    
    return ilr, xlr, ℓ
end
