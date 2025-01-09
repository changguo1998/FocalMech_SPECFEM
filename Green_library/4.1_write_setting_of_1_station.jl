using DelimitedFiles, SeisTools, Printf, TOML

function proc_multi(r1::AbstractFloat, r2::AbstractFloat, rh::AbstractFloat)
    l = Tuple{Int,Int,Int,Int}[]
    for p1 = 1:8, p2 = 1:8, m1 = 1:10, m2 = 1:10
        h1 = r1/(p1*m1*8)
        h2 = r2/(p2*m2*8)
        if (max(h1, h2)/min(h1, h2) < 1.2) && (abs(sqrt(h1*h2)-rh) < 0.5)
            push!(l, (p1, m1, p2, m2))
        end
    end
    np = map(_x->_x[1]*_x[3], l)
    return l[argmax(np)]
end

function locatepoint(x::Real, xs::AbstractVector{<:Real}, ERRoutofrange::Bool=true)
    if ERRoutofrange && ((x < minimum(xs)) || (x > maximum(xs)))
        error(string(x)*" Out of Model Range")
    end
    i = findfirst(>(x), xs)
    if isnothing(i)
        i = length(xs)
        h = 1.0
    elseif i == 1
        i = 2
        h = 0.0
    else
        h = (x - xs[i-1]) / (xs[i] - xs[i-1])
    end
    return (i, h)
end

function bilinear(x, i, j, p, q)
    return x[i-1, j-1] * (1.0 - p) * (1.0 - q) +
           x[i-1, j] * (1.0 - p) * q +
           x[i, j-1] * p * (1.0 - q) +
           x[i, j] * p * q
end

cutmodeldir = "/public/home/guochang/environment/chuandian_cutmodel"

(topolat, topolon, topo) = let
    # io = open(joinpath(@__DIR__, "topography.bin"))
    io = open(joinpath(cutmodeldir, "topography.bin"))
    nlon = read(io, Int32)
    nlat = read(io, Int32)
    tlon = zeros(Float32, nlon)
    tlat = zeros(Float32, nlat)
    read!(io, tlon)
    read!(io, tlat)
    dat = zeros(Float32, nlon, nlat)
    read!(io, dat)
    close(io)
    (tlat, tlon, dat)
end;

# station_range = readdlm("../../dat/statistic/station_area.csv"; comments=true)
station_range = let
    buf = filter(!startswith("#"), readlines("../../dat/statistic/station_area.csv"))
    tl = map(buf) do l
        m = split(l, ' ', keepempty=false)
        (String(m[1]), parse.(Float64, m[2:end]))
    end
    tab = Matrix{Any}(undef, length(buf), 7)
    for i = eachindex(tl)
        tab[i, 1] = tl[i][1]
        tab[i, 2:7] = tl[i][2][1:6]
    end
    tab
end

margin = 10.0
pml = [20.0, 20.0, 20.0]
dt = 0.1
dx = 0.5
dz = 0.1
fddh = 0.5
# dx = 5.0
# dz = 1.0
npts = 1800
rt = 1.0
zmin = -8.0
zmax = 100.0
garea_dx = 0.5
gdx = 2.0
gdz = 1.0
gzmax = 50.0
mesh_h = 2.0


C0 = pi*6371.0/180.0

# basicRootDir = "/public/home/guochang/semsimu/chuandian2/dat/glib"
basicRootDir = "../../dat/glib"
# Threads.@threads for s = axes(station_range, 1)
# for s = [57, 80]
let
    s = parse(Int, ARGS[1])
    (tlat, t_p) = locatepoint(station_range[s, 2], topolat)
    (tlon, t_q) = locatepoint(station_range[s, 3], topolon)
    dp = bilinear(topo, tlon, tlat, t_q, t_p)

    minx = (station_range[s, 4] - station_range[s, 2])*C0 - margin |> floor
    maxx = (station_range[s, 5] - station_range[s, 2])*C0 + margin |> ceil
    miny = (station_range[s, 6] - station_range[s, 3])*C0*cosd(station_range[s, 2]) - margin |> floor
    maxy = (station_range[s, 7] - station_range[s, 3])*C0*cosd(station_range[s, 2]) + margin |> ceil
    modelrange = [maxx-minx, maxy-miny, zmax]
    t = proc_multi(modelrange[1], modelrange[2], mesh_h)

    config = Dict("wkdir"          => abspath(basicRootDir, station_range[s, 1]),
                "modelrange"     => modelrange,
                "pml"            => pml,
                "srcloc"         => [-minx, -miny, ceil(dp/dz)*dz],
                "srclat"         => station_range[s, 2],
                "srclon"         => station_range[s, 3],
                "srcdep"         => ceil(dp/dz)*dz,
                "dt"             => dt,
                "npts"           => npts,
                "amp"            => 1e10,
                "risetime"       => rt,
                "modelfile"      => abspath(basicRootDir, station_range[s, 1], "model.bin"),
                "nproc_xi"       => t[3],
                "nproc_eta"      => t[1],
                "xi_multiply"    => t[4],
                "eta_multiply"   => t[2],
                "step_multiply"  => 20,
                "irregular_mesh" => true,
                "surface"        => [8],
                "tlib_h"         => fddh,
                "receiver"       => Dict{String,Any}[])

    narea_x = round(Int, (station_range[s, 5] - station_range[s, 4])/garea_dx)
    narea_y = round(Int, (station_range[s, 7] - station_range[s, 6])/garea_dx)
    p = 0
    for i = 1:narea_x, j = 1:narea_y
        p += 1
        lat1 = station_range[s, 4] + garea_dx * (i-1)
        lat2 = station_range[s, 4] + garea_dx * i
        lon1 = station_range[s, 6] + garea_dx * (j-1)
        lon2 = station_range[s, 6] + garea_dx * j
        x1 = (lat1 - station_range[s, 2])*C0 - minx
        x2 = (lat2 - station_range[s, 2])*C0 - minx
        y1 = (lon1 - station_range[s, 3])*C0*cosd(station_range[s, 2]) - miny
        y2 = (lon2 - station_range[s, 3])*C0*cosd(station_range[s, 2]) - miny
        push!(config["receiver"], Dict{String, Any}(
            "id" => p,
            "n" => [x1, gdx, x2],
            "e" => [y1, gdx, y2],
            "d" => [zmin, gdz, gzmax]
        ))
    end

    mkpath("../../dat/glib/"*station_range[s, 1])
    open(abspath("../../dat/glib/", station_range[s, 1], "setting.toml"), "w") do io
        TOML.print(io, config; sorted=true)
    end

    cmd1 = Cmd(["julia", "-t", "63", joinpath(cutmodeldir, "genmodel.jl"),
        "-D", string(zmax),
        "-H", string(-zmin),
        "-M", "2",
        "-O", abspath("../../dat/glib/", station_range[s, 1], "model.bin"),
        "-R", join([minx, maxx, miny, maxy], '/'),
        "-N", join(round.(Int, [modelrange[1]/dx, modelrange[2]/dx, (zmax-zmin)/dz]), '/'),
        "-S", string(station_range[s, 2])*"/"*string(station_range[s, 3])])
    println(cmd1)
    try
        if !isdir(abspath(abspath("../../dat/glib/", station_range[s, 1], "model.bin")))
            run(cmd1)
        end
    catch err
        println(err)
    end
end