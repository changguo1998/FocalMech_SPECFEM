using SeisTools

statable = let
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

# Nsta = size(statable, 1)

# dist_north = zeros(Nsta)
# dist_east = zeros(Nsta)
# for i = axes(statable, 1)
#     dist_north[i] = SeisTools.Geodesy.distance(statable[i, 4], 0.0, 
#         statable[i, 5], 0.0)
#     dist_east[i] = SeisTools.Geodesy.distance(statable[i, 4], statable[i, 6], 
#         statable[i, 4], statable[i, 7])
# end

# Nbatch = 5
# NeventperBatch = ceil(Int, Nsta/Nbatch)

# azmap = zeros(Nsta, Nsta)
# distmap = zeros(Nsta, Nsta)
# for i = axes(statable, 1), j = axes(statable, 1)
#     azmap[i, j] = SeisTools.Geodesy.azimuth(statable[i, 2], statable[i, 3],
#         statable[j, 2], statable[j, 3])
#     distmap[i, j] = SeisTools.Geodesy.distance(statable[i, 2], statable[i, 3],
#         statable[j, 2], statable[j, 3])*1e-3
# end

# # seedlist = [122, 31, 139, 126, 133]
# seedlist = [98, 104,80, 86, 145]

# batchid = zeros(Int, Nsta)
# for i = 1:Nbatch
#     iseed = seedlist[i]
#     idxs = sortperm(distmap[iseed, :])
#     npick = 0
#     for j = idxs
#         if iszero(batchid[j])
#             batchid[j] = i
#             npick += 1
#         end
#         if npick == NeventperBatch || !any(iszero, batchid)
#             break
#         end
#     end
# end

# for ib = 1:Nbatch
#     open("batch_$(ib)_note.md", "w") do io
#         println(io, '|', join(["hash", "S", "T", "D", "E", "N", "G", "C", "K"], '|'), '|')
#         println(io, "|:-:"^9, '|')
#         for i = findall(==(i), batchid)
#             println(io, "|", statable[i, 1], "|"^9)
#         end
#     end
# end

# for ib = 1:Nbatch
#     open("batch_$(ib).txt", "w") do io
#         # println(io, '|', join(["hash", "S", "T", "D", "E", "N", "G", "C", "K"], '|'), '|')
#         # println(io, "|:-:"^9, '|')
#         for i = findall(==(ib), batchid)
#             # println(io, "|", statable[i, 1], "|"^9)
#             println(io, statable[i, 1])
#         end
#     end
# end

batchfiles = filter(startswith("batch"), readdir())
stationlist = String[]
batchids = Int[]
for ib in eachindex(batchfiles)
    l = readlines(batchfiles[ib])
    append!(stationlist, l)
    append!(batchids, fill(ib, length(l)))
end

spm = map(statable[:, 1]) do s
    findfirst(==(s), stationlist)
end

stationlist = stationlist[spm]
batchid = batchids[spm]

clist = [:red, :blue, :green, :gray70, :purple]

fig = Figure(resolution=(800,1000))
ax = Axis(fig[1, 1];xticks=98:0.5:107, yticks=22:0.5:33)
vlines!(collect(98:0.5:107), color=:gray80, width=1)
hlines!(collect(22:0.5:33), color=:gray80, width=1)
for i = 1:maximum(batchid)
    fidx = batchid .== i
    scatter!(Float64.(statable[fidx, 3]), Float64.(statable[fidx, 2]),
        color=clist[i], marker=:utriangle, markersize=16)
end

text!(Float64.(statable[:, 3]), Float64.(statable[:, 2]);
        text=String.(statable[:, 1]), align=(:center, :top))
