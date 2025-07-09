using Printf, Downloads

function _srtm_get_par(io::IO, T::Type)
    t = readline(io)
    l = split(t, ' ', keepempty=false)
    return parse(T, l[2])
end

function _read_asc_header(fn::String)
    io = open(fn)
    ncol = _srtm_get_par(io, Int)
    nrow = _srtm_get_par(io, Int)
    olon = _srtm_get_par(io, Float64)
    olat = _srtm_get_par(io, Float64)
    dd = _srtm_get_par(io, Float64)
    nanvalue = _srtm_get_par(io, Float64)
    close(io)
    return (nrow, ncol, olat, olon, dd, nanvalue)
end

function readasc(fn::AbstractString)
    tls = readlines(fn)
    ncol = parse(Int, tls[1][15:end])
    nrow = parse(Int, tls[2][15:end])
    cx = parse(Float64, tls[3][15:end])
    cy = parse(Float64, tls[4][15:end])
    cs = parse(Float64, tls[5][15:end])
    nanvalue = parse(Float64, tls[6][15:end])
    dat = zeros(Float64, nrow, ncol)
    for ir = 1:nrow
        l = tls[ir+6]
        sl = split(l)
        for ic = 1:ncol
            dat[ir, ic] = parse(Float64, sl[ic])
            if dat[ir, ic] == nanvalue
                dat[ir, ic] = NaN
            end
        end
    end
    return (cx, cy, cs, dat)
end

#
ascfilelist = readdir(joinpath(@__DIR__, "unpack"))
ascheaders = map(_f->_read_asc_header(joinpath(@__DIR__, "unpack", _f)), ascfilelist)

lonidx = map(ascfilelist) do f
    return parse(Int, f[6:7])
end

latidx = map(ascfilelist) do f
    return parse(Int, f[9:10])
end

latlist = sort(unique(latidx))
lonlist = sort(unique(lonidx))
minlat = minimum(getindex.(ascheaders, 3))
maxlat = maximum(getindex.(ascheaders, 3).+getindex.(ascheaders, 1) .* getindex.(ascheaders, 5))
minlon = minimum(getindex.(ascheaders, 4))
maxlon = maximum(getindex.(ascheaders, 4).+getindex.(ascheaders, 2) .* getindex.(ascheaders, 5))

# Modify the 6 variables below according to downloaded block
nblock_lat = length(latlist)
iblock_lat_start = minimum(latlist)
lat_range = (minlat, maxlat)
nblock_lon = length(lonlist)
iblock_lon_start = minimum(lonlist)
lon_range = (minlon, maxlon)

buffer = fill(Float32(10000.0), nblock_lat*5999+1, nblock_lon*5999+1)

for ilon = 1:nblock_lon, ilat = 1:nblock_lat
    fname = @sprintf("srtm_%d_%02d.asc", ilon+iblock_lon_start-1, ilat+iblock_lat_start-1)
    println(fname)
    local fpath = joinpath(@__DIR__, "unpack", fname)
    if !isfile(fpath)
        continue
    end
    slon = (ilon-1)*5999
    slat = (ilat-1)*5999
    (_, _, _, dat) = readasc(joinpath(@__DIR__, "unpack", fname))
    buffer[slat.+axes(dat, 1), slon.+axes(dat, 2)] .= dat
end

data = -permutedims(buffer)./1000.0
reverse!(data, dims=2)
lat = range(start=lat_range[1], stop=lat_range[2], length=size(data, 2))
lon = range(start=lon_range[1], stop=lon_range[2], length=size(data, 1))

open("topography.bin", "w") do io
    write(io, Int32(length(lon)))
    write(io, Int32(length(lat)))
    write(io, Float32.(lon))
    write(io, Float32.(lat))
    write(io, Float32.(data))
end
