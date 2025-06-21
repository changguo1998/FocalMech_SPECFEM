using DelimitedFiles

(table,_) = readdlm("SWChinaCVM-V2.0-main/SWChinaCVMv2.0.txt.wrst.sea_level", header=true)

lons = unique(table[:, 1]) |> sort
lats = unique(table[:, 2]) |> sort
deps = unique(table[:, 3]) |> sort

minlon = lons[1]
minlat = lats[1]
mindep = deps[1]

dlon = minimum(diff(lons))
maxlonid = round(Int, (lons[end]-minlon)/dlon) + 1
lonid = zeros(Int, maxlonid)
for i = eachindex(lons)
    lid = round(Int, (lons[i]-minlon)/dlon)+1
    lonid[lid] = i
end

dlat = minimum(diff(lats))
maxlatid = round(Int, (lats[end]-minlat)/dlat) + 1
latid = zeros(Int, maxlatid)
for i = eachindex(lats)
    lid = round(Int, (lats[i]-minlat)/dlat)+1
    latid[lid] = i
end

dd = minimum(diff(deps))
maxdepid = round(Int, (deps[end]-mindep)/dd) + 1
depid = zeros(Int, maxdepid)
for i = eachindex(deps)
    did = round(Int, (deps[i]-mindep)/dd) + 1
    depid[did] = i
end

vp = zeros(Float32, length(deps), length(lons), length(lats))
vs = zeros(Float32, length(deps), length(lons), length(lats))

for i = axes(table, 1)
    itx = round(Int, (table[i, 2]-minlat)/dlat)+1
    ity = round(Int, (table[i, 1]-minlon)/dlon)+1
    itz = round(Int, (table[i, 3]-mindep)/dd)+1
    ix = latid[itx]
    iy = lonid[ity]
    iz = depid[itz]
    vp[iz, iy, ix] = table[i, 4]
    vs[iz, iy, ix] = table[i, 5]
end

open("CVM2.bin", "w") do io
    write(io, Int32(length(lats)))
    write(io, Int32(length(lons)))
    write(io, Int32(length(deps)))
    write(io, Float32.(lats))
    write(io, Float32.(lons))
    write(io, Float32.(deps))
    write(io, vp)
    write(io, vs)
end
