using NCDatasets

"""
    carpenter2020buoyancy(nav_file, ctd_file)
Create a buoyancy field from NetCDF from Carpenter et al 2020
"""
function carpenter2020buoyancy(nav_file, ctd_file)
    ds = NCDataset(nav_file)
    time_pos = ds["time_pos"][:]
    east_pos = ds["east_pos"][:]
    north_pos = ds["north_pos"][:]
    close(ds)
    
    ds = NCDataset(ctd_file)
    time_ctd = ds["time_ctd"][:]
    depth_ctd = ds["depth_ctd"][:]
    ρ_ctd = ds["potential_density_ctd"][:]
    close(ds)
    
    along_track_distance = similar(east_pos)
    dx = diff(east_pos)
    dy = diff(north_pos)
    ds = @. sqrt(dx^2 + dy^2)
    along_track_distance[2:end] = cumsum(ds)
    ρ_ref = 1005.35
    
    pos_to_ctd_index = Iterators.map(t->findfirst(time_ctd .>= t), time_pos)
    
    x_ctd = similar(time_pos)
    z_ctd = similar(time_pos)
    b_ctd = similar(time_pos)

    println("Reading $(length(pos_to_ctd_index)) positions")
    for (pos_index, ctd_index) in enumerate(pos_to_ctd_index)
        rem(pos_index, 100) == 1 && print("$pos_index, $ctd_index\r")
        x_ctd[pos_index] = along_track_distance[pos_index]
        z_ctd[pos_index] = -depth_ctd[ctd_index]
        b_ctd[pos_index] = (ρ_ref-ρ_ctd[ctd_index]) / ρ_ref * 9.81
    end
    println()
    Nx, Nz = 512, 128
    
    xs_face = range(minimum(x_ctd), maximum(x_ctd), Nx+1)
    zs_face = range(minimum(z_ctd) - 5, 0, Nz+1)
    xs_center = (xs_face[2:end] .+ xs_face[1:end-1]) ./ 2
    zs_center = (zs_face[2:end] .+ zs_face[1:end-1]) ./ 2
    
    Δx = xs_center[2] - xs_center[1]
    Δz = zs_center[2] - zs_center[1]
    
    bs = zeros(eltype(b_ctd), Nx, Nz)
    ws = zeros(eltype(b_ctd), Nx, Nz)
    maximum_depth = zeros(eltype(b_ctd), Nx)
    println("Interpolating $(length(x_ctd)) readings")
    for i_ctd in eachindex(x_ctd)
        rem(i_ctd, 100) == 1 && print("$i_ctd\r")
        x, z, b = x_ctd[i_ctd], z_ctd[i_ctd], b_ctd[i_ctd]
        i = convert(typeof(i_ctd), 1 + (x - xs_center[1]) ÷ Δx)
        j = convert(typeof(i_ctd), 1 + (z - zs_center[1]) ÷ Δz)
        maximum_depth[i] = max(maximum_depth[i], -z)
        for Δi in -4:4, Δj in -4:4
            if (1 <= i+Δi <= Nx) && (1 <= j+Δj <= Nz)
                δx = (x - xs_center[i+Δi]) / Δx
                δz = (z - zs_center[j+Δj]) / Δz
                w = exp(-(δx^2 + δz^2) / 2)
                
                bs[i+Δi, j+Δj] += w * b
                ws[i+Δi, j+Δj] += w
            end
        end
    end

    maximum_depth = map(eachindex(maximum_depth)) do i
        maximum(maximum_depth[max(i-4, 1):min(i+4, Nx)])
    end
    
    println()
    
    ε = exp(-8)
    bs = @. bs / (ws .+ ε)
    
    println("Creating field")
    grid = RectilinearGrid(; 
        x = (xs_face[1], xs_face[end]), 
        z = (zs_face[1], zs_face[end]), 
        size = (Nx, Nz),
        topology = (Bounded, Flat, Bounded)
    )
    
    carpenter_buoyancy = Field{Center, Center, Center}(grid)
    set!(carpenter_buoyancy, bs)
    fill_halo_regions!(carpenter_buoyancy)

    H = Field{Center, Center, Nothing}(grid)
    set!(H, maximum_depth)
    fill_halo_regions!(H)
    
    return carpenter_buoyancy, H
end

@inline function carpenter2020preprocessfunc(i, j, k, grid, carpenter_buoyancy, maximum_depth, x_carpenter_center, b_min)
    from_loc = (Center(), nothing, Center())
    from_grid = carpenter_buoyancy.grid
    
    x, z = node(i, j, k, grid, from_loc...)
    x_carpenter = x_carpenter_center - x
    z_carpenter = z
    to_node = (x_carpenter, z_carpenter)

    z_min = -Oceananigans.interpolate(x_carpenter, maximum_depth, (Center(), nothing, nothing), from_grid)
    z_carpenter < z_min && return b_min
    return Oceananigans.interpolate(to_node, carpenter_buoyancy, from_loc, from_grid)
end

"""
    carpenter2020preprocess(carpenter_buoyancy, maximum_depth)
Process field `carpenter_buoyancy` to match our simulation output.
"""
function carpenter2020preprocess(carpenter_buoyancy, maximum_depth)
    carpenter_grid = carpenter_buoyancy.grid
    carpenter_x = xnodes(carpenter_buoyancy)
    
    # Region of interest, note reflection
    x_right = 39e3
    x_left = 51e3
    x_center = (x_right + x_left) / 2

    # Parameters...
    Lh = abs(x_right - x_left)
    Lx = Lh
    Lz = -znode(1, carpenter_buoyancy.grid, Face())
    Nx = Integer(Lx ÷ (carpenter_x[2] - carpenter_x[1]))
    Nh = Nx ÷ 2
    Ny = 1
    Nz = size(carpenter_buoyancy.grid)[3]
    Ly = Lx / Nx
    sp = (; Nx, Ny, Nz, Lx, Ly, Lz, Nh, Lh)
    
    # New grid, 3D so it works with any badly written things I made...
    grid = RectilinearGrid(;
        x = (-Lx/2, Lx/2),
        y = (-Ly/2, Ly/2),
        z = (-Lz, 0),
        size = (Nx, Ny, Nz),
        topology = (Bounded, Periodic, Bounded)
    )
    b_min = minimum(carpenter_buoyancy)
    b = Field(KernelFunctionOperation{Center, Nothing, Center}(carpenter2020preprocessfunc, grid, carpenter_buoyancy, maximum_depth, x_center, b_min))
    compute!(b)
    fill_halo_regions!(b)

    b_surface = Field(ZSlice(b, 0.0))
    b_min = minimum(b_surface)
    b_max = maximum(b_surface)
    sp = (sp..., Δb = b_max - b_min)
    
    return b, b_surface, sp
end
