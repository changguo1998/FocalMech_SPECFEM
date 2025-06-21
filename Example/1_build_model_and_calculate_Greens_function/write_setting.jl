#
# station_range format:
# station_name,station_lat,station_lon,glib_min_lat,glib_max_lat,glib_min_lon,glib_max_lon
#
range_north = 50.0
range_east = 50.0
range_depth = 30.0
pml_thickness = 5.0

config = Dict(
    "wkdir" => abspath(pwd()),
    "modelrange" => [range_north, range_east, range_depth],
    "pml" => [pml_thickness, pml_thickness, pml_thickness],
    "srcloc" => [-minx, -miny, ceil(dp / dz) * dz],
    "srclat" => station_range[s, 2],
    "srclon" => station_range[s, 3],
    "srcdep" => ceil(dp / dz) * dz,
    "dt" => dt,
    "npts" => npts,
    "amp" => 1e10,
    "risetime" => rt,
    "modelfile" => abspath(basicRootDir, station_range[s, 1], "model.bin"),
    "nproc_xi" => t[3],
    "nproc_eta" => t[1],
    "xi_multiply" => t[4],
    "eta_multiply" => t[2],
    "step_multiply" => 20,
    "irregular_mesh" => true,
    "surface" => [8],
    "tlib_h" => fddh,
    "receiver" => Dict{String,Any}[]
)

narea_x = round(Int, (station_range[s, 5] - station_range[s, 4]) / garea_dx)
narea_y = round(Int, (station_range[s, 7] - station_range[s, 6]) / garea_dx)
p = 0
for i = 1:narea_x, j = 1:narea_y
    p += 1
    lat1 = station_range[s, 4] + garea_dx * (i - 1)
    lat2 = station_range[s, 4] + garea_dx * i
    lon1 = station_range[s, 6] + garea_dx * (j - 1)
    lon2 = station_range[s, 6] + garea_dx * j
    x1 = (lat1 - station_range[s, 2]) * C0 - minx
    x2 = (lat2 - station_range[s, 2]) * C0 - minx
    y1 = (lon1 - station_range[s, 3]) * C0 * cosd(station_range[s, 2]) - miny
    y2 = (lon2 - station_range[s, 3]) * C0 * cosd(station_range[s, 2]) - miny
    push!(config["receiver"], Dict{String,Any}(
        "id" => p,
        "n" => [x1, gdx, x2],
        "e" => [y1, gdx, y2],
        "d" => [zmin, gdz, gzmax]
    ))
end

mkpath("../../dat/glib/" * station_range[s, 1])
open(abspath("../../dat/glib/", station_range[s, 1], "setting.toml"), "w") do io
    TOML.print(io, config; sorted=true)
end

cmd1 = Cmd(["julia", "-t", "63", joinpath(cutmodeldir, "genmodel.jl"),
    "-D", string(zmax),
    "-H", string(-zmin),
    "-M", "2",
    "-O", abspath("../../dat/glib/", station_range[s, 1], "model.bin"),
    "-R", join([minx, maxx, miny, maxy], '/'),
    "-N", join(round.(Int, [modelrange[1] / dx, modelrange[2] / dx, (zmax - zmin) / dz]), '/'),
    "-S", string(station_range[s, 2]) * "/" * string(station_range[s, 3])])
println(cmd1)
try
    if !isdir(abspath(abspath("../../dat/glib/", station_range[s, 1], "model.bin")))
        run(cmd1)
    end
catch err
    println(err)
end
