using CairoMakie, Mmap

fpath = abspath(@__DIR__, "topography.bin")
@info "Loading data from $fpath..."

io = open(fpath)
nlon = read(io, Int32)
nlat  = read(io, Int32)
lons = zeros(Float32, nlon)
lats = zeros(Float32, nlat)
read!(io, lons)
read!(io, lats)
data = Mmap.mmap(io, Matrix{Float32}, (nlon, nlat), 4 * (2 + nlon + nlat))

r0 = 1000 / max(nlon, nlat)

resampled_data = zeros(round(Int, nlon * r0 * 2), round(Int, nlat * r0 * 2))

indexs_along_lon = round.(Int, range(1, nlon, size(resampled_data, 1)))
indexs_along_lat = round.(Int, range(1, nlat, size(resampled_data, 2)))
for j = eachindex(indexs_along_lat), i = eachindex(indexs_along_lon)
    resampled_data[i, j] = -data[indexs_along_lon[i], indexs_along_lat[j]]
end

close(io)

# for test
# resampled_data[1:round(Int, size(resampled_data, 1)/2), 1:round(Int, size(resampled_data, 1)/2)] .= 10.0

fig = Figure(size=(round(Int, nlon * r0), round(Int, nlat * r0)));
ax  = Axis(fig[1, 1], title="Topography", xlabel="Longitude (°E)", ylabel="Latitude (°N)");

hp = heatmap!(ax, lons[indexs_along_lon], lats[indexs_along_lat], resampled_data;
    colormap=:topo, highclip=:red, colorrange=(-10.0, 9.9));
Colorbar(fig[1,2], hp);
save("topography.png", fig);
