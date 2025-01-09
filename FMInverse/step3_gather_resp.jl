include("env.jl")

respdir = joinpath(PTH["root"], "dat", "resp")
if !isdir(respdir)
    mkpath(respdir)
end

for batchdir = (PTH["batch1"], PTH["batch2"])
    for edir in readdir(batchdir)
        Threads.@threads for f in readdir(joinpath(batchdir, edir))
            if !startswith(f, "RESP")
                continue
            end
            if isfile(joinpath(respdir, f*"_"*edir))
                continue
            end
            cp(joinpath(batchdir, edir, f), joinpath(respdir, f*"_"*edir))
        end
    end
end

respraw = readdir(respdir)
channellist = map(rname->split(rname, '_')[1], respraw) |> unique

resp_uniq = joinpath(PTH["root"], "dat", "RESP")
mkpath(resp_uniq)

for cha = channellist
    local fs = filter(startswith(cha), respraw)
    local sameid = zeros(Int, length(fs))
    for i = eachindex(fs), j = eachindex(fs)
        if i > j
            continue
        end
        if sameid[j] > 0
            continue
        end
        local t1 = read(joinpath(respdir, fs[i]), String)
        local t2 = read(joinpath(respdir, fs[j]), String)
        if t1 == t2
            sameid[j] = i
        end
    end
    local fus = unique(sameid)
    if length(fus) == 1
        if !isfile(joinpath(resp_uniq, cha))
            cp(joinpath(respdir, fs[fus[1]]), joinpath(resp_uniq, cha))
        end
    else
        for iu = fus
            if isfile(joinpath(resp_uniq, cha*"_$iu"))
                continue
            end
            cp(joinpath(respdir, fs[iu]), joinpath(resp_uniq, cha*"_$iu"))
        end
    end
end

buffer = reduce(vcat, map(readlines, readdir(resp_uniq, join=true)))
open(io->foreach(l->println(io, l), buffer), joinpath(PTH["root"], "dat", "RESP.all"), "w")
