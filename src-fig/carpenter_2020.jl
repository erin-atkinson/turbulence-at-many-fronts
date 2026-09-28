using NCDatasets

ds = NCDataset("NAV_data.nc")
time_pos = ds["time_pos"][:]
east_pos = ds["east_pos"][:]
north_pos = ds["north_pos"][:]
close(ds)

ds = NCDataset("CTD_data.nc")
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

(x_ctd, z_ctd, b_ctd) = begin
    pos_to_ctd_index = Iterators.map(t->findfirst(time_ctd .>= t), time_pos)
    
    x_ctd = similar(time_pos)
    z_ctd = similar(time_pos)
    b_ctd = similar(time_pos)

    for (pos_index, ctd_index) in enumerate(pos_to_ctd_index)
        print("$pos_index, $ctd_index\r")
        x_ctd[pos_index] = along_track_distance[pos_index]
        z_ctd[pos_index] = -depth_ctd[ctd_index]
        b_ctd[pos_index] = (σ_max-σ_ctd[ctd_index]) / σ_max * 9.81
    end
    
    x_ctd, z_ctd, b_ctd
end

Nx, Nz = 1024, 128

xs = range(minimum(x_ctd), maximum(x_ctd), Nx)
zs = range(minimum(z_ctd), maximum(z_ctd), Nz)

Δx = xs[2] - xs[1]
Δz = zs[2] - zs[1]
bs = zeros(eltype(b_ctd), Nx, Nz)
ws = zeros(eltype(b_ctd), Nx, Nz)

for i_ctd in eachindex(x_ctd)
    x, z, b = x_ctd[i_ctd], z_ctd[i_ctd], b_ctd[i_ctd]
    i = convert(typeof(i_ctd), 1 + (x - xs[1]) ÷ Δx)
    j = convert(typeof(i_ctd), 1 + (z - zs[1]) ÷ Δz)
    for Δi in -5:5, Δj in -5:5
        if (1 < i+Δi < Nx) && (1 < j+Δj < Nz)
            δx = (x - xs[i+Δi]) / Δx
            δz = (z - zs[j+Δj]) / Δz
            w = exp(-(δx^2 + δz^2) / 2)
            
            bs[i+Δi, j+Δj] += w * b
            ws[i+Δi, j+Δj] += w
        end
    end
end
bs = @. bs / (ws .+ 1e-8)

b_max = maximum(bs)
b_min = b_max - 0.002
colorrange = (b_min, b_max)
colormap = reverse(to_colormap(:deep))

fig = Figure()
ax = Axis(fig[1, 1];
    xlabel = L"\text{Along-track distance} / \text{km}",
    ylabel = L"z / \text{m}",
    limits = (38, 52, -75, 0)
)
contourf!(xs ./ 1000, zs, bs; colormap, levels=range(b_min, b_max, 10), extendlow=:black)
fig