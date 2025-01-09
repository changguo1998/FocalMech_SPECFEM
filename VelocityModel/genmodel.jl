using Pkg
Pkg.activate(abspath(@__DIR__, ".."))
using ArgumentProcessor

include("semmodel.jl")

const r0 = π * 6371.0 / 180.0

addopt!("output"; abbr="O", fmt=" %s", default=" " * joinpath(pwd(), "model.bin"),
        help="Output model file")
addopt!("station"; abbr="S", fmt=" %f/%f", required=true, help="lat/lon")
# addopt!("radius"; abbr="R", fmt=" %f", required=true, help="radius of model, unit in kilometer(km)")
addopt!("range"; abbr="R", fmt=" %f/%f/%f/%f", required=true, help="horizontal range Nmin/Nmax/Emin/Emax")
addopt!("maxelevation"; abbr="H", fmt=" %f", default=" 9.0", help="elevation of top surface of model")
addopt!("depth"; abbr="D", fmt=" %f", required=true, help="depth of bottom surface of model")
addopt!("sample"; abbr="N", fmt=" %d/%d/%d", required=true,
        help="n sample in each direction. nlat/nlon/ndep or nx/ny/nz")
addopt!("model"; abbr="M", fmt=" %d", required=true, help="1 CVM1.0; 2 CVM2.0")

input = ArgumentProcessor.parse(ARGS)

# radius = input.radius
hrange = input.range
modelfile = abspath(input.output)

stalat = input.station[1]
stalon = input.station[2]

# minlat = stalat - radius/r0
# maxlat = stalat + radius/r0
minlat = stalat + hrange[1]/r0
maxlat = stalat + hrange[2]/r0
avelat = (minlat + maxlat) / 2.0
nlat = input.sample[1]

# minlon = stalon - radius/r0/cosd(avelat)
# maxlon = stalon + radius/r0/cosd(avelat)
minlon = stalon + hrange[3]/r0/cosd(avelat)
maxlon = stalon + hrange[4]/r0/cosd(avelat)
nlon = input.sample[2]

mindep = -input.maxelevation
maxdep = input.depth
ndep = input.sample[3]

println("lat range: ", minlat, " ", maxlat)
println("lon range: ", minlon, " ", maxlon)

# if minlat < 21.0
#     error("Out of Range, minlat too little")
# end

# if maxlat > 34.0
#     error("Out of Range, maxlat too large")
# end

# if minlon < 97.0
#     error("Out of Range, minlon too little")
# end

# if maxlon > 108.0
#     error("Out of Range, maxlon too large")
# end

function writemodelbin(fn::AbstractString, vp::AbstractArray, vs::AbstractArray, rho::AbstractArray,
                       topz::Real, dx::Real, dy::Real, dz::Real, airmask::AbstractArray{Bool}=Bool[])
    (v0, dv, nsp, materialhash, materialtable) = discretematerial((vp, vs, rho), (200, 200, 50), .!airmask)
    materialgrid = field2material((vp, vs, rho), v0, dv, materialhash, airmask)

    open(fn, "w") do io
        # * field
        write(io, Int32(length(v0)))
        write(io, Float32.(v0))
        write(io, Float32.(dv))
        write(io, Int32.(nsp))
        write(io, Int32(size(materialtable, 2)))
        write(io, Float32.(materialtable))
        # * mesh
        write(io, Int32(size(vp, 3)))
        write(io, Int32(size(vp, 2)))
        write(io, Int32(size(vp, 1)))
        write(io, Float32(dx))
        write(io, Float32(dy))
        write(io, Float32(dz))
        write(io, Float32(topz))
        # * model
        write(io, Int32.(materialgrid))
    end
    return nothing
end

(latlist, lonlist, deplist, mvp, mvs) = let
    if input.model == 1
        io = open(joinpath(@__DIR__, "SWChinaCVM1.0.sea_level.bin"))
    elseif input.model == 2
        io = open(joinpath(@__DIR__, "SWChinaCVM2.0.sea_level.bin"))
    else
        error("model type not correct. see help for more information.")
    end
    nlat = read(io, Int32)
    nlon = read(io, Int32)
    ndep = read(io, Int32)
    lats = zeros(Float32, nlat)
    lons = zeros(Float32, nlon)
    deps = zeros(Float32, ndep)
    read!(io, lats)
    read!(io, lons)
    read!(io, deps)
    vp = zeros(Float32, ndep, nlon, nlat)
    vs = zeros(Float32, ndep, nlon, nlat)
    read!(io, vp)
    read!(io, vs)
    close(io)
    (lats, lons, deps, vp, vs)
end

(topolat, topolon, topo) = let
    io = open(joinpath(@__DIR__, "topography.bin"))
    nlon = read(io, Int32)
    nlat = read(io, Int32)
    tlon = zeros(Float32, nlon)
    tlat = zeros(Float32, nlat)
    read!(io, tlon)
    read!(io, tlat)
    dat = zeros(Float32, nlon, nlat)
    read!(io, dat)
    (tlat, tlon, dat)
end

if (minlat < minimum(latlist)) || (minlat < minimum(topolat))
    error("Out of Range, minlat too little")
end

if (maxlat > maximum(latlist)) || (maxlat > maximum(topolat))
    error("Out of Range, maxlat too large")
end

if (minlon < minimum(lonlist)) || (minlon < minimum(topolon))
    error("Out of Range, minlon too little")
end

if (maxlon > maximum(lonlist)) || (maxlon > maximum(topolon))
    error("Out of Range, maxlon too large")
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

function triplelinear(x, i, j, k, p, q, r)
    if any(iszero, x[i-1:i, j-1:j, k-1:k])
        error("Velocity is Zero")
    end
    return x[i-1, j-1, k-1] * (1.0 - p) * (1.0 - q) * (1.0 - r) +
           x[i-1, j-1, k] * (1.0 - p) * (1.0 - q) * r +
           x[i-1, j, k-1] * (1.0 - p) * q * (1.0 - r) +
           x[i-1, j, k] * (1.0 - p) * q * r +
           x[i, j-1, k-1] * p * (1.0 - q) * (1.0 - r) +
           x[i, j-1, k] * p * (1.0 - q) * r +
           x[i, j, k-1] * p * q * (1.0 - r) +
           x[i, j, k] * p * q * r
end

vp = zeros(Float32, ndep, nlon, nlat)
vs = zeros(Float32, ndep, nlon, nlat)
rho = fill(Float32(2.7), ndep, nlon, nlat)
isair = falses(ndep, nlon, nlat)
Threads.@threads for idx in CartesianIndices(vp)
    (idep, ilon, ilat) = idx.I
    lat = (ilat - 1) / (nlat - 1) * (maxlat - minlat) + minlat
    lon = (ilon - 1) / (nlon - 1) * (maxlon - minlon) + minlon
    dep = (idep - 1) / (ndep - 1) * (maxdep - mindep) + mindep
    (tlat, t_p) = locatepoint(lat, topolat)
    (tlon, t_q) = locatepoint(lon, topolon)
    ptopo = bilinear(topo, tlon, tlat, t_q, t_p)

    isair[idx] = dep < ptopo
    if !isair[idx]
        (mlat, p) = locatepoint(lat, latlist)
        (mlon, q) = locatepoint(lon, lonlist)
        (mdep, r) = locatepoint(dep, deplist, false)
        vp[idx] = triplelinear(mvp, mdep, mlon, mlat, r, q, p)
        vs[idx] = triplelinear(mvs, mdep, mlon, mlat, r, q, p)
    end
end

if input.model == 1
    println("Using velocity model: CVM1.0")
else
    println("Using velocity model: CVM2.0")
end
println("Model write to file: ", modelfile)
writemodelbin(modelfile, vp, vs, rho, mindep,
              r0 * (maxlat - minlat) / (nlat - 1),
              r0 * cosd(avelat) * (maxlon - minlon) / (nlon - 1),
              (maxdep - mindep) / (ndep - 1), isair)
