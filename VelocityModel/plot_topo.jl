using CairoMakie

fpath = abspath(@__DIR__, "topography.bin")
@info "Loading data from $fpath..."

io = open(fpath)
nlon = read(io, Int32)
nlat  = read(io, Int32)
lons = zeros(Float32, nlon)
lats = zeros(Float32, nlat)
data = zeros(Float32, nlon, nlat)

read!(io, lons)
read!(io, lats)
read!(io, data)
close(io)

data = -data

r0 = 1000 / max(nlon, nlat)

fig = Figure(size=(round(Int, nlon * r0), round(Int, nlat * r0)));
ax  = Axis(fig[1, 1], title="Topography", xlabel="Longitude (°E)", ylabel="Latitude (°N)");

hp = heatmap!(ax, lons, lats, data; colormap=:topo, colorrange=(-8.0, 8.0));
Colorbar(fig[1,2], hp);
save("topography.png", fig);
