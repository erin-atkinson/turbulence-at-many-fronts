using MAT

"""
    mckie_2024_buoyancy(filename; ρ_name="dens", p_name="P", lat_name="lat", lon_name="lon", flip_x=false)
Create a buoyancy field from MAT data from McKie 2024
"""
function mckie_2024_buoyancy(filename; ρ_name="dens", p_name="P", lat_name="lat", lon_name="lon", flip_x=false)
    b = matopen(filename) do file
        println(keys(file))
        ρ = reverse(transpose(read(file, ρ_name)); dims=2)
        @. 9.81 * (1020 - ρ) / 1020
    end;
    z = matopen(filename) do file
        p = read(file, p_name)[:]
        reverse(-p)
    end;
    lat = matopen(filename) do file
        read(file, lat_name)[:]
    end;
    lon = matopen(filename) do file
        read(file, lon_name)[:]
    end;
    
    R = 6_371_000
    dlat = [lat[i] - lat[max(i-1, 1)] for i in eachindex(lat)] .* (π / 180)
    dlon = [lon[i] - lon[max(i-1, 1)] for i in eachindex(lon)] .* (π / 180)
    dx = @. R * sqrt(cos(lat)^2 * dlon^2 + dlat^2)
    x = cumsum(dx)
    if flip_x
        x = maximum(x) .- reverse(x)
        reverse!(b; dims=1)
    end
    b_fixed = fix_NaN(b)
    
    Lz = maximum(abs, z)
    xs_center = x
    Δxᶜ = diff(xs_center)
    
    xs_face = zeros(eltype(xs_center), length(xs_center)+1)
    xs_face[1] = xs_center[1] - (Δxᶜ[1]/2)
    xs_face[end] = xs_center[end] + (Δxᶜ[end]/2)
    xs_face[2:end-1] = xs_center[1:end-1] .+ (Δxᶜ/2)
    
    grid = RectilinearGrid(; 
        x = xs_face, 
        z = (-Lz, 0.0), 
        size = (length(xs_center), length(z)),
        topology = (Bounded, Flat, Bounded)
    )
    
    b_field = CenterField(grid)
    set!(b_field, b_fixed)
    fill_halo_regions!(b_field)
    
    return b_field
end

function fix_NaN(b)
    b_fixed = similar(b)
    
    for i in axes(b, 1)
        k1 = findfirst(!isnan, b[i, :]) + 1
        k2 = findlast(!isnan, b[i, :]) - 5
        for k in axes(b, 2)
            b_fixed[i, k] = b[i, min(k2, max(k1, k))]
        end
    end
    
    return filt(b_fixed, 2, 1)
end

@inline function mckie_2024_preprocess_func(i, j, k, grid, mckie_buoyancy, x_center)
    from_loc = (Center(), nothing, Center())
    from_grid = mckie_buoyancy.grid
    
    x, z = node(i, j, k, grid, from_loc...)
    x_mckie = x + x_center
    z_mckie = z
    to_node = (x_mckie, z_mckie)

    return Oceananigans.interpolate(to_node, mckie_buoyancy, from_loc, from_grid)
end

"""
    mckie_2024_preprocess(mckie_buoyancy, x_right, x_left)
Process field `mckie_buoyancy` to match our simulation output.
"""
function mckie_2024_preprocess(mckie_buoyancy, x_left, x_right)
    mckie_grid = mckie_buoyancy.grid
    mckie_x = xnodes(mckie_buoyancy)
    
    # Region of interest
    x_center = (x_right + x_left) / 2

    # Parameters...
    Lh = abs(x_right - x_left)
    Lx = Lh
    Lz = -znode(1, mckie_buoyancy.grid, Face())
    Nx = Integer(Lx ÷ (mckie_x[2] - mckie_x[1]))
    Nh = Nx ÷ 2
    Ny = 1
    Nz = size(mckie_buoyancy.grid)[3]
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
    b_min = minimum(mckie_buoyancy)
    b = Field(KernelFunctionOperation{Center, Nothing, Center}(mckie_2024_preprocess_func, grid, mckie_buoyancy, x_center))
    compute!(b)
    fill_halo_regions!(b)

    b_surface = Field(ZSlice(b, 0.0))
    b_min = minimum(b_surface)
    b_max = maximum(b_surface)
    sp = (sp..., Δb = b_max - b_min, f=3.76e-5)
    
    return b, b_surface, sp
end
