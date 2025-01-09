using Printf, SeisTools, DelimitedFiles

function checkavaiable(meshtype::Matrix{Int}, i1::Int, i2::Int, j1::Int, j2::Int)
    all(>(0), meshtype[i1:i2,j1:j2])
end

function mesh2range(meshtype::Matrix{Int})
    L1 = size(meshtype, 1)
    L2 = size(meshtype, 2)
    Ti = maximum(meshtype, dims=2) |> vec
    Tj = maximum(meshtype, dims=1) |> vec
    I1 = findfirst(==(2), Ti)
    I2 = findlast(==(2), Ti)
    J1 = findfirst(==(2), Tj)
    J2 = findlast(==(2), Tj)
    area = 0.0
    irange = [0,0,0,0]
    for i1 = I1:I2-1, j1 = J1:J2-1
        for i2 = i1+1:I2, j2 = j1+1:J2
            if !checkavaiable(meshtype, i1, i2, j1, j2)
                continue
            end
            S = (i2-i1+1)*(j2-j1+1)
            if sum(meshtype[i1:i2,j1:j2]) < 2 * S
                continue
            end
            if area < S
                area = S
                irange[1] = i1
                irange[2] = i2
                irange[3] = j1
                irange[4] = j2
            end
        end
    end
    return irange
end

C0 = pi*6371.0/180.0

radius = 1.8*C0
margin = 0.4*C0
step = 0.5
table = let
    txts = readlines("../../dat/statistic/stationlist_cb.txt")
    t = Matrix{Any}(undef, (length(txts), 3))
    for i = eachindex(txts)
        l = split(txts[i])
        t[i, 1] = l[1]
        t[i, 2] = parse(Float64, l[2])
        t[i, 3] = parse(Float64, l[3])
    end
    t
end

(cvm2, _) = readdlm("../../dat/cvm2/SWChinaCVMv2.0.txt.wrst.sea_level"; header=true)
latlist = sort(unique(cvm2[:, 2]))
lonlist = sort(unique(cvm2[:, 1]))
deplist = sort(unique(cvm2[:, 3]))
depcount = zeros(Int, length(latlist), length(lonlist))
for i = axes(cvm2, 1)
    ilon = findfirst(abs.(lonlist .- cvm2[i, 1]) .< 1e-3)
    ilat = findfirst(abs.(latlist .- cvm2[i, 2]) .< 1e-3)
    idep = findfirst(abs.(deplist .- cvm2[i, 3]) .< 1e-3)
    depcount[ilat, ilon] += 1
end
filledmap = depcount .== length(deplist)

lats = 22.0:step:33.0
lons = 98.0:step:107.0
meshtype = zeros(Int, length(lats), length(lons))
meshmask = trues(length(lats), length(lons))
for i = axes(meshtype, 1), j = axes(meshtype, 2)
    for p = axes(filledmap, 1), q = axes(filledmap, 2)
        if SeisTools.Geodesy.distance(lats[i], lons[j], latlist[p], lonlist[q])/1000 < margin
            meshmask[i, j] &= filledmap[p, q]
        end
    end
end

meshRange = zeros(size(table, 1), 4)

for s = axes(table, 1)
    for i = axes(meshtype, 1), j = axes(meshtype, 2)
        # flag = max(abs(lons[j] - table[s, 3])/180.0*pi*6371.0*cosd(table[s, 2]),
        #     abs(lats[i] - table[s, 2])/180.0*pi*6371.0) >= radius
        flag = radius <= max(
            SeisTools.Geodesy.distance(table[s, 2], table[s, 3], table[s, 2], lons[j]),
            SeisTools.Geodesy.distance(table[s, 2], table[s, 3], lats[i], table[s, 3])
        )*0.001
        if meshmask[i, j]
            if flag
                meshtype[i, j] = 1
            else
                meshtype[i, j] = 2
            end
        else
            meshtype[i, j] = 0
        end
    end
    r = mesh2range(meshtype)
    meshRange[s, 1] = lats[r[1]]
    meshRange[s, 2] = lats[r[2]]
    meshRange[s, 3] = lons[r[3]]
    meshRange[s, 4] = lons[r[4]]
end

idxs = Int[]
for i = axes(table, 1)
    if table[i, 2] > meshRange[i, 1] &&
        table[i, 2] < meshRange[i, 2] &&
        table[i, 3] > meshRange[i, 3] &&
        table[i, 3] < meshRange[i, 4]
        push!(idxs, i)
    end
end

open("../../dat/statistic/station_area.csv", "w") do io
    println(io, "#hash lat lon minlat maxlat minlon maxlon")
    for i = idxs
        @printf(io, "%s %7.4f %8.4f %4.1f %4.1f %5.1f %5.1f\n", table[i,:]..., meshRange[i,:]...)
    end
end
