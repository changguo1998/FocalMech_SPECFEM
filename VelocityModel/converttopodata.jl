using Printf

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

buffer = zeros(Float32, 3*6000-2, 3*6000-2)

for ilon = 1:3, ilat = 1:3
    fname = @sprintf("srtm_%d_%02d.asc", ilon+55, ilat+5)
    println(fname)
    slon = (ilon-1)*5999
    slat = (ilat-1)*5999
    (_, _, _, dat) = readasc(joinpath(@__DIR__, "unpack", fname))
    buffer[slat.+axes(dat, 1), slon.+axes(dat, 2)] .= dat
end

data = -permutedims(buffer)./1000.0
reverse!(data, dims=2)
lat = range(start=20.0, stop=35.0, length=size(data, 2))
lon = range(start=95.0, stop=110.0, length=size(data, 1))

open("topography.bin", "w") do io
    write(io, Int32(length(lon)))
    write(io, Int32(length(lat)))
    write(io, Float32.(lon))
    write(io, Float32.(lat))
    write(io, Float32.(data))
end
