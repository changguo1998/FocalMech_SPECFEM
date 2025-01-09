using Printf, DelimitedFiles, SeisTools, SHA

GD = SeisTools.Geodesy
# = ========================================================
tab0 = [
    let 
        t = readdlm("../../dat/statistic/stationlist_v1.txt")
        map(i->(t[i,1], round(t[i, 2]; digits=4), round(t[i, 3]; digits=4)), axes(t, 1))
    end;
    let
        t = readdlm("../../dat/statistic/stationlist_unpack.txt")
        map(i->(t[i,1], round(t[i, 2]; digits=4), round(t[i, 3]; digits=4)), axes(t, 1))
    end
]

# = ========================================================
tab1 = eltype(tab0)[]
for _t = tab0
    if !(_t in tab1)
        push!(tab1, _t)
    end
end

# = =======================================================
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

tab2 = filter(tab1) do t
    ilat = findlast(latlist .< t[2])
    ilon = findlast(lonlist .< t[3])
    !isnothing(ilat) && !isnothing(ilon) &&
    latlist[end] > t[2] && lonlist[end] > t[3] &&
    filledmap[ilat,ilon] && filledmap[ilat+1,ilon] &&
    filledmap[ilat,ilon+1] && filledmap[ilat+1,ilon+1]
end

# = ========================================================
tab3 = sort(tab2, by=x->String(x[1][4:end]))

# = ========================================================
dmap = zeros(size(tab3, 1), size(tab3, 1))
for i = axes(dmap, 1), j = axes(dmap, 2)
    dmap[i, j] = GD.distance(tab3[i][2], tab3[i][3], tab3[j][2], tab3[j][3])/1e3
end

idxlist = CartesianIndex[]
for idx = findall(<=(20.0), dmap)
    if (idx.I[1] >= idx.I[2]) || (tab3[idx.I[1]][1][4:end] != tab3[idx.I[2]][1][4:end])
        continue
    end
    push!(idxlist, idx)
end
flags = trues(length(tab3))
for idx = idxlist
    if startswith(tab3[idx.I[1]][1], "XX")
        flags[idx.I[1]] = false
    end
    if startswith(tab3[idx.I[2]][1], "XX")
        flags[idx.I[2]] = false
    end
end
tab4 = tab3[flags]

# = ========================================================
dmap = zeros(size(tab4, 1), size(tab4, 1))
for i = axes(dmap, 1), j = axes(dmap, 2)
    dmap[i, j] = GD.distance(tab4[i][2], tab4[i][3], tab4[j][2], tab4[j][3])/1e3
end

idxlist = CartesianIndex[]
for idx = findall(<=(5.0), dmap)
    if (idx.I[1] >= idx.I[2])
        continue
    end
    push!(idxlist, idx)
end
flags = trues(length(tab4))
lowprec(x) = abs(x - round(x*100.0)*0.01) < 1e-5
for idx = idxlist
    if lowprec(tab4[idx.I[1]][2]) && lowprec(tab4[idx.I[1]][3])
        flags[idx.I[1]] = false
    end
    if lowprec(tab4[idx.I[2]][2]) && lowprec(tab4[idx.I[2]][3])
        flags[idx.I[2]] = false
    end
end
tab5 = tab4[flags]

# = ========================================================
dmap = zeros(size(tab5, 1), size(tab5, 1))
for i = axes(dmap, 1), j = axes(dmap, 2)
    dmap[i, j] = GD.distance(tab5[i][2], tab5[i][3], tab5[j][2], tab5[j][3])/1e3
end

l1 = Int[]
cblist = filter(i->i.I[1]<i.I[2], findall(<(5.0), dmap))
for _i = cblist
    append!(l1, _i.I)
end
sort!(l1)
l2 = filter(i->!(i in l1), eachindex(tab5))

hashcoor(lat,lon) = reinterpret(UInt8, [Float64(lat), Float64(lon)]) |> 
    sha256 |> x->string.(x[1:4]; base=16, pad=2) |> join |> uppercase

tab6 = Tuple{String,Float64,Float64,Vector{String}}[]

for i = l2
    lat = tab5[i][2]
    lon = tab5[i][3]
    push!(tab6, (hashcoor(lat, lon), lat, lon, String[]))
end
for cb = cblist
    i, j = cb.I
    lat = round(tab5[i][2]; digits=4)
    lon = round(tab5[i][3]; digits=4)
    push!(tab6, (hashcoor(lat, lon), lat, lon, String[]))
end

for i = eachindex(tab5)
    for j = eachindex(tab6)
        if GD.distance(tab5[i][2], tab5[i][3], tab6[j][2], tab6[j][3]) <= 5000.0
            if !(tab5[i][1] in tab6[j][4])
                push!(tab6[j][4], tab5[i][1])
            end
        end
    end
end

sort!(tab6; by=x->String(x[4][1][4:end]))

open("../../dat/statistic/stationlist_cb.txt", "w") do io
    for t in tab6
        @printf(io, "%s %7.4f %8.4f", t[1], t[2], t[3])
        for l in t[4]
            print(io, " ", l)
        end
        print(io, "\n")
    end
end

#=
fig = Figure(resolution=(800, 1000))
ax = Axis(fig[1, 1]; autolimitaspect=1)
text!(map(t->t[3:-1:2], tab5); text=getindex.(tab5, 1))
scatter!(map(t->t[3:-1:2], tab5), marker=:utriangle, markersize=16)
linelist = Tuple{Float64,Float64,String}[]
for idx = findall(<=(5.0), dmap)
    if (idx.I[1] >= idx.I[2])
        continue
    end
    push!(linelist, (
        (tab5[idx.I[1]][3]+tab5[idx.I[2]][3])/2, 
        (tab5[idx.I[1]][2]+tab5[idx.I[2]][2])/2,
        @sprintf("%.4f", dmap[idx])
    ))
    scatterlines!([tab5[idx.I[1]][3], tab5[idx.I[2]][3]], [tab5[idx.I[1]][2], tab5[idx.I[2]][2]], 
        color=:red, linewidth=2, markersize=6)
end
text!(map(t->t[1:2], linelist); text=getindex.(linelist, 3))
=#
