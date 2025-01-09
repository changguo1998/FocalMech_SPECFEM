using SeisTools, DelimitedFiles
include("env.jl")
include("lib/io.jl")
include("lib/lib.jl")

SAC = SeisTools.SAC
Geo = SeisTools.Geodesy
DP = SeisTools.DataProcess

function findnearest(lat::Real, lon::Real, slist::Vector)
    dists = map(s->Geo.distance(lat, lon, s[2], s[3])*1e-3, slist)
    return findmin(dists)
end

function networkstr(s::AbstractString)
    t = split(s, '.')
    return String(t[1])
end

function hasunmatch(d1::Dict, d2::Dict, kl::Vector{String})
    for k in kl
        if typeof(d1[k]) <: Real
            if abs(d1[k]-d2[k]) > 1e-3
                return true
            end
        else
            if d1[k] != d2[k]
                return true
            end
        end
    end
    return false
end

function waveinfo(hdr::Dict)
    reftime = SAC.DateTime(hdr)
    bt = reftime + Millisecond(round(Int, hdr["b"]*1e3))
    dt = Millisecond(round(Int, hdr["delta"]*1e3))
    et = bt + hdr["npts"] * dt
    return (bt, dt, et)
end

function is_same_event(edir1::AbstractString, edir2::AbstractString)
    sacfs1 = readdir(edir1, join=true)
    sacfs2 = readdir(edir2, join=true)
    nmatch = 0
    spairs = Tuple{String,String}[]
    for s1 in sacfs1, s2 in sacfs2
        push!(spairs, (String(s1), String(s2)))
    end
    flags = falses(length(spairs))
    Threads.@threads for ip in eachindex(spairs)
        (s1, s2) = spairs[ip]
        h1 = open(SAC.readhead, s1)
        h2 = open(SAC.readhead, s2)
        if hasunmatch(h1, h2, ["stla", "stlo", "delta"])
            continue
        end
        (bt1, dt1, et1) = waveinfo(h1)
        (bt2, dt2, et2) = waveinfo(h2)
        if (bt1 > et2) || (bt2 > et1)
            continue
        end
        d1 = open(SAC.read, s1)
        d2 = open(SAC.read, s2)
        bt = max(bt1, bt2)
        et = min(et1, et2)
        cbt = bt + (et - bt) / 2
        w1 = DP.cut(d1.data, bt1, cbt, Minute(1), dt1)
        w2 = DP.cut(d2.data, bt2, cbt, Minute(1), dt2)
        res = sum(abs, w1[2]-w2[2])
        if res > 1.0e-3
            continue
        end
        flags[ip] = true
    end
    return count(flags) > 0
end

function dirname2dt(e::AbstractString)
    if startswith(e, "2")
        estr = e
    else
        estr = String(e[4:15])
    end
    return parseymd(estr)
end

# ==========================================================
source = PTH["dat_filted"];
target = PTH["dat_selectstation"];

# # load
pr1 = PhaseIO.readphasereport("processed_phasereport.txt");
pr2 = let
    (t, _) = readdlm(joinpath(PTH["root"], "dat/catalog_need_waveform.txt"), header=true)
    map(axes(t, 1)) do i
        PhaseIO.Event(
            (region=("XX", ""),
            par=(DateTime(t[i, 1], "y/m/d,H:M:S"), t[i, 3], t[i, 2], t[i, 4], t[i, 5]),
            phase=PhaseIO.Phase[])
        )
    end
end;
pr = [pr1; pr2];
etime = map(_p->_p.par[1], pr);

stationlist = map(readlines(PTH["stationlist"])) do l
    b = split(l, ' '; keepempty=false)
    hsh = b[1]
    lat = parse(Float64, b[2])
    lon = parse(Float64, b[3])
    nms = b[4:end]
    (hsh, lat, lon, nms)
end;

@info "index"
event_copy_schedule = Vector{String}[];
event_candindates = readdir(source);
flag_classified = falses(length(event_candindates));
for i = eachindex(event_candindates)
    if flag_classified[i]
        continue
    end
    local evts = String[]
    push!(evts, event_candindates[i])
    flag_classified[i] = true
    for j = eachindex(event_candindates)
        global event_candindates, flag_classified
        if i >= j
            continue
        end
        if abs(dirname2dt(event_candindates[i]) - dirname2dt(event_candindates[j])) > Minute(5)
            continue
        end
        if flag_classified[j]
            continue
        end
        if is_same_event(joinpath(source, event_candindates[i]), joinpath(source, event_candindates[j]))
            push!(evts, String(event_candindates[j]))
            flag_classified[j] = true
        end
    end
    push!(event_copy_schedule, evts)
    # println("($(count(flag_classified[i:end])))", i, "/", length(event_candindates))
end

open("duplicated_event.txt", "w") do io
    for es in event_copy_schedule
        for e in es
            print(io, " ", e)
        end
        println(io, "")
    end
end

@info "copy"
for elist = event_copy_schedule
    local etgtsacdir = joinpath(target, elist[1], "sac")
    mkpath(etgtsacdir)
    local stacplist = Dict{String,String}()
    for e in elist
        for f = readdir(joinpath(source, e))
            local td, id
            local sacsrcpath = joinpath(source, e, f)
            local hdr = open(SAC.readhead, sacsrcpath)
            (td, id) = findnearest(hdr["stla"], hdr["stlo"], stationlist)
            if td > 2.0
                continue
            end
            local tag = join([networkstr(stationlist[id][4][1]), hdr["kstnm"], hdr["kcmpnm"]], ".")
            if tag in keys(stacplist)
                continue
            end
            stacplist[tag] = sacsrcpath
        end
    end
    for k in keys(stacplist)
        local hdr, dat, td, id
        (hdr, dat) = open(SAC.read, stacplist[k])
        (td, id) = findnearest(hdr["stla"], hdr["stlo"], stationlist)
        hdr["knetwk"] = networkstr(stationlist[id][4][1])
        hdr["khole"] = "00"
        fnewname = SAC.standardname(hdr; standard="cenc")
        fnewpath = joinpath(etgtsacdir, fnewname)
        SAC.write(fnewpath, hdr, dat)
    end
end

@info "write phase report"
for e = readdir(target)
    et = dirname2dt(e)
    ie = argmin(abs.(etime .- et))
    if abs(etime[ie]-et) < Minute(1)
        if (pr[ie].par[2] > 22.5) &&
            (pr[ie].par[2] < 32.7) &&
            (pr[ie].par[3] > 98.5) &&
            (pr[ie].par[3] < 105.5) &&
            ( -0.8 * (pr[ie].par[3]-101.0) + 22.5 < pr[ie].par[2])
            PhaseIO.writephasereport(joinpath(target, e, "phase.txt"), [pr[ie]])
        end
    end
end

@info "delete empty"
for e = readdir(target)
    if isempty(readdir(joinpath(target, e, "sac"))) ||
        (!isfile(joinpath(target, e, "phase.txt")))
        rm(joinpath(target, e), recursive=true)
    end
end
