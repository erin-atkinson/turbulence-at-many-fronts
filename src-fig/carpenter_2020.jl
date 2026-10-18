using NCDatasets

function carpenter2020buoyancy(nav_file, ctd_file)
    ds = NCDataset(nav_file)
    time_pos = ds["time_pos"][:]
    east_pos = ds["east_pos"][:]
    north_pos = ds["north_pos"][:]
    close(ds)
    
    ds = NCDataset(ctd_file)
    time_ctd = ds["time_ctd"][:]
    depth_ctd = ds["depth_ctd"][:]
    σ_ctd = ds["potential_density_ctd"][:]
    close(ds)
    
    along_track_distance = similar(east_pos)
    dx = diff(east_pos)
    dy = diff(north_pos)
    ds = @. sqrt(dx^2 + dy^2)
    along_track_distance[2:end] = cumsum(ds)
    σ_max = maximum(σ_ctd)
    
    pos_to_ctd_index = Iterators.map(t->findfirst(time_ctd .>= t), time_pos)
    
    x_ctd = similar(time_pos)
    z_ctd = similar(time_pos)
    b_ctd = similar(time_pos)

    println("Reading $(length(pos_to_ctd_index)) positions")
    for (pos_index, ctd_index) in enumerate(pos_to_ctd_index)
        rem(pos_index, 100) == 1 && print("$pos_index, $ctd_index\r")
        x_ctd[pos_index] = along_track_distance[pos_index]
        z_ctd[pos_index] = -depth_ctd[ctd_index]
        b_ctd[pos_index] = (σ_max-σ_ctd[ctd_index]) / σ_max * 9.81
    end
    println()
    
    Nx, Nz = 512, 64
    
    xs_face = range(minimum(x_ctd), maximum(x_ctd), Nx+1)
    zs_face = range(minimum(z_ctd), 0, Nz+1)
    
    xs_center = (xs_face[2:end] .+ xs_face[1:end-1]) ./ 2
    zs_center = (zs_face[2:end] .+ zs_face[1:end-1]) ./ 2
    
    Δx = xs_center[2] - xs_center[1]
    Δz = zs_center[2] - zs_center[1]
    
    bs = zeros(eltype(b_ctd), Nx, Nz)
    ws = zeros(eltype(b_ctd), Nx, Nz)
    println("Interpolating $(length(x_ctd)) readings")
    for i_ctd in eachindex(x_ctd)
        rem(i_ctd, 100) == 1 && print("$i_ctd\r")
        x, z, b = x_ctd[i_ctd], z_ctd[i_ctd], b_ctd[i_ctd]
        i = convert(typeof(i_ctd), 1 + (x - xs_center[1]) ÷ Δx)
        j = convert(typeof(i_ctd), 1 + (z - zs_center[1]) ÷ Δz)
        for Δi in -10:10, Δj in -5:5
            if (1 <= i+Δi <= Nx) && (1 <= j+Δj <= Nz)
                δx = (x - xs_center[i+Δi]) / Δx
                δz = (z - zs_center[j+Δj]) / Δz
                w = exp(-(δx^2 + δz^2) / 2)
                
                bs[i+Δi, j+Δj] += w * b
                ws[i+Δi, j+Δj] += w
            end
        end
    end
    println()
    
    # Frontal region between 44 and 50 km
    bs = @. bs / (ws .+ 1e-8)
    
    is = filter(i -> 44_000 < xs_center[i] < 50_000, eachindex(xs_center))
    i_right = argmin(i -> bs[i, Nz], is)
    i_left = argmax(i -> bs[i, Nz], is)
    
    x₀ = xs_center[i_right]
    Δb = abs(bs[i_right] - bs[i_left])
    b_max = maximum(bs[i_left:i_right, :])
    bs = bs .- b_max
    
    println("Creating field")
    # On a 3D grid just so any naive functions I made work ok
    grid = RectilinearGrid(; 
        x = xs_face, 
        y = [-Δx/2, Δx/2], 
        z = zs_face, 
        size = (Nx, 1, Nz),
        topology = (Bounded, Periodic, Bounded)
    )
    
    field = Field{Center, Nothing, Center}(grid)
    set!(field, reverse(bs; dims=1)) # match our simulations
    fill_halo_regions!(field)

    # Just from paper and by eye
    H = 60
    Lz = 75
    f = 1.2e-4
    L = sqrt(Δb * H) / f
    
    sp = (; Δb, f, H, L, Lz)
    
    return field, sp
end