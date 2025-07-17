#
# station_range format:
# station_name,station_lat,station_lon,glib_min_lat,glib_max_lat,glib_min_lon,glib_max_lon
#

using TOML

function proc_multi(r1::AbstractFloat, r2::AbstractFloat, rh::AbstractFloat, ncores::Integer)
    l = Tuple{Int,Int,Int,Int}[]
    m_max = ceil(Int, max(r1, r2)/rh/8)
    for p1 = 1:ncores, p2 = 1:ncores
        if p1*p2 > ncores
            continue
        end
        for m1 = 1:m_max, m2 = 1:m_max
            h1 = r1/(p1*m1*8)
            h2 = r2/(p2*m2*8)
            if (max(h1, h2)/min(h1, h2) < 1.2) && (abs(sqrt(h1*h2)-rh) < 0.5)
                push!(l, (p1, m1, p2, m2))
            end
        end
    end
    np = map(_x->_x[1]*_x[3], l)
    return l[argmax(np)]
end

C0 = pi*6371.0/180.0

cutmodeldir = abspath(@__DIR__, "..", "..", "VelocityModel")
velmodelpath = abspath(@__DIR__, "..", "..", "VelocityModel", "CVM2.bin")

pml_thickness = 20.0
station_range = ["STATION" 30.0 100.0 29.0 31.0 99.0 101.0]

margin = 10.0
dt = 0.1
fddh = 0.5
dx = 5.0
dz = 0.1
npts = 900
rt = 4.0
zmin = -8.0
zmax = 100.0
garea_dx = 0.5
gdx = 2.0
gdz = 1.0
gzmax = 30.0
mesh_h = 5.0

s = 1 # index of station_range, useful when there are multiple stations.

minx = (station_range[s, 4] - station_range[s, 2])*C0 - margin |> floor
maxx = (station_range[s, 5] - station_range[s, 2])*C0 + margin |> ceil
miny = (station_range[s, 6] - station_range[s, 3])*C0*cosd(station_range[s, 2]) - margin |> floor
maxy = (station_range[s, 7] - station_range[s, 3])*C0*cosd(station_range[s, 2]) + margin |> ceil
modelrange = [maxx-minx, maxy-miny, zmax]
t = proc_multi(modelrange[1], modelrange[2], mesh_h, round(Int, Sys.CPU_THREADS/2))

wkdir = abspath(pwd(), station_range[s, 1])

config = Dict(
    "wkdir" => wkdir,
    "modelrange" => modelrange,
    "pml" => [pml_thickness, pml_thickness, pml_thickness],
    "srcloc" => [-minx, -miny, 0.0],
    "srclat" => station_range[s, 2],
    "srclon" => station_range[s, 3],
    "srcdep" => 0.0,
    "dt" => dt,
    "npts" => npts,
    "amp" => 1e10,
    "risetime" => rt,
    "modelfile" => abspath(wkdir, "model.bin"),
    "nproc_xi" => t[3],
    "nproc_eta" => t[1],
    "xi_multiply" => t[4],
    "eta_multiply" => t[2],
    "step_multiply" => 10,
    "irregular_mesh" => true,
    "surface" => [8],
    "tlib_h" => fddh,
    "receiver" => Dict{String,Any}[]
)

narea_x = round(Int, (station_range[s, 5] - station_range[s, 4]) / garea_dx)
narea_y = round(Int, (station_range[s, 7] - station_range[s, 6]) / garea_dx)
p = 0
for i = 1:narea_x, j = 1:narea_y
    global p += 1
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

mkpath(wkdir)
open(abspath(wkdir, "setting.toml"), "w") do io
    TOML.print(io, config; sorted=true)
end

cmd1 = Cmd(["julia", "-t", string(Sys.CPU_THREADS-1), joinpath(cutmodeldir, "genmodel.jl"),
    "-D", string(zmax),
    "-H", string(-zmin),
    "-M", velmodelpath,
    "-O", abspath(wkdir, "model.bin"),
    "-R", join([minx, maxx, miny, maxy], '/'),
    "-N", join(round.(Int, [modelrange[1] / dx, modelrange[2] / dx, (zmax - zmin) / dz]), '/'),
    "-S", string(station_range[s, 2]) * "/" * string(station_range[s, 3])])
println(cmd1)
try
    if !isdir(abspath(abspath(wkdir, "model.bin")))
        run(cmd1)
    end
catch err
    println(err)
end
