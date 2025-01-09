using Dates, SeisTools, DelimitedFiles

SAC = SeisTools.SAC
Geo = SeisTools.Geodesy

include("env.jl")
include("lib/io.jl")
include("lib/lib.jl")

target = PTH["dat_filted"]
pr = PhaseIO.readphasereport("processed_phasereport.txt")
ecatalog = let
    (t, _) = readdlm(joinpath(PTH["root"], "dat/catalog_need_waveform.txt"), header=true)
    t
end
etime = [map(_p->_p.par[1], pr); DateTime.(ecatalog[:, 1], "y/m/d,H:M:S")]
elons = [map(_p->_p.par[3], pr); Float64.(ecatalog[:, 2])]
elats = [map(_p->_p.par[2], pr); Float64.(ecatalog[:, 3])]

event_not_match = String[]
for source = (PTH["batch1"], PTH["batch2"])
    for e in readdir(source)
        if startswith(e, "2")
            estr = e
        else
            estr = String(e[4:15])
        end
        ie = argmin(abs.(etime .- parseymd(estr)))
        if abs(etime[ie]-parseymd(estr)) > Minute(1)
            push!(event_not_match, joinpath(source, e)*":"*string(etime[ie]))
            continue
        end
        elat = elats[ie]
        elon = elons[ie]
        mkpath(joinpath(target, e))
        for f in readdir(joinpath(source, e))
            if endswith(f, "SAC")
                hdr = open(SAC.readhead, joinpath(source, e, f))
                slat = hdr["stla"]
                slon = hdr["stlo"]
                dist = Geo.distance(elat, elon, slat, slon)*1e-3
                if dist < 200.0
                    if !isfile(joinpath(target, e, f))
                        cp(joinpath(source, e, f), joinpath(target, e, f); force=true)
                    end
                end
            end
        end
    end
end

open(io->foreach(l->println(io, l), event_not_match), "rerange_not_match.txt", "w")

# open("rerange_not_match.txt", "w") do io
#     for l in event_not_match
#         println(io, leading_zeros)
#     end
# end