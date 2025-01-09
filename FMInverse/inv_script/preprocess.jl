using DelimitedFiles, SeisTools, TOML, Dates, JuliaSourceMechanism,
    SourceMechUI, Printf, JLD2, Statistics, LinearAlgebra

Geo = SeisTools.Geodesy

include("../env.jl")

include(joinpath(PTH["invsrc"], "lib/io.jl"))
include(joinpath(PTH["invsrc"], "lib/lib.jl"))
include(joinpath(@__DIR__, "multistage_lib.jl"))
include("/home/guochang/Projects/NLLOCshell/NLLOCshell.jl")
include("/home/guochang/Projects/glibSEMshell/include/io.jl")

misfitmodules = [XCorr, Polarity]

# eventname = "201304200802479"
# eventname = "202106140402371"

eventname = ARGS[1]

# dataroot = abspath("..")
# println("dataroot: ", dataroot)

# eventname = splitpath(dataroot)[end]
# dataroot = "/data/guochang/chuandian2/dat/selectstation/201304200802479"
dataroot = abspath(PTH["dat_selectstation"], eventname)


glibindex = let
    t = readlines(joinpath(PTH["glib"], "index.txt"))
    map(t) do l
        s = split(l)
        (s[1], parse(Float64, s[2]), parse(Float64, s[3]))
    end
end
glats = getindex.(glibindex, 2)
glons = getindex.(glibindex, 3)
pr=PhaseIO.readphasereport(joinpath(dataroot, "phase.txt"))
ereport=pr[1]
event = Dict("longitude" => ereport.par[3],
        "latitude" => ereport.par[2],
        "depth" => ereport.par[4] > 30.0 ? 30.0 : ereport.par[4],
        "magnitude" => ereport.par[5],
        "origintime" => ereport.par[1] - Hour(8))

algorithm = Dict("misfit" => [m.tags[1] for m in misfitmodules],
        "searchdepth" => event["depth"],
        "weight" => ones(length(misfitmodules)))

stations = buildstationconfiguration(dataroot, event)

TIME_DEFAULT = DateTime(2000)
Eareacode = area_code(event["latitude"], event["longitude"])
# = station infomation

C0 = pi*6371.0/180.0

glibexistflag = falses(length(stations))
for i = eachindex(glibexistflag)
    local s = stations[i]
    # sid = map(glibindex) do _s
    #     SeisTools.Geodesy.distance(s["meta_lat"], s["meta_lon"], _s[2], _s[3])
    # end |> argmin
    # sdist = SeisTools.Geodesy.distance(s["meta_lat"], s["meta_lon"], glibindex[sid][2], glibindex[sid][3])
    (sid, sdist) = nearest_station_in_glib(s["meta_lat"], s["meta_lon"], glats, glons)
    if sdist > 1000.0
        glibexistflag[i] = false
        continue
    end
    glibpath = joinpath(PTH["glib"], glibindex[sid][1], Eareacode*".bin")
    if isfile(glibpath)
        x = (event["latitude"] - s["meta_lat"]) * C0
        y = (event["longitude"] - s["meta_lon"]) * C0 * cosd(s["meta_lat"])
        io = open(glibpath)
        (nx, ny, nz, x0, y0, z0, dx, dy, dz, _, _, _, _,
            _, _, _, _, _, _, _, _, _, _, _) = cglib_readhead(io)
        close(io)
        if (x > (x0+(nx-1)*dx)) || (x < x0) ||
            (y > (y0+(ny-1)*dy)) || (y < y0)
            glibexistflag[i] = false
        else
            glibexistflag[i] = true
        end
    else
        glibexistflag[i] = false
    end
end

stations = stations[glibexistflag]

for s in stations
    # println(s["station"])
    sacf = SeisTools.SAC.read(normpath(dataroot, "sac", s["meta_file"]))
    # - trim window
    reft = SeisTools.SAC.DateTime(sacf.hdr)
    bt = event["origintime"] - Minute(5)
    et = bt + Minute(15)
    # et = reft + msecond(sacf.hdr["delta"] * sacf.hdr["npts"])
    s["base_trim"] = [bt, et]
    ps = filter(ereport.phase) do p
        (p.network == s["network"]) &&
        (p.station == s["station"]) &&
        (p.type[1] == 'P')
    end
    if isempty(ps)
        push!(s["phases"], Dict(
            "type" => "P",
            "at" => TIME_DEFAULT
        ))
    else
        at_test = minimum(getfield.(ps, :time))
        if at_test > Time(ereport.par[1])
            _atp = DateTime(Date(ereport.par[1]), at_test) - Hour(8)
        else
            _atp = DateTime(Date(ereport.par[1])+Day(1), at_test) - Hour(8)
        end
        push!(s["phases"], Dict(
            "type" => "P",
            "at" => _atp
        ))
    end
    ps = filter(ereport.phase) do p
        (p.network == s["network"]) &&
        (p.station == s["station"]) &&
        (p.type[1] == 'S')
    end
    if isempty(ps)
        push!(s["phases"], Dict(
            "type" => "S",
            "at" => TIME_DEFAULT
        ))
    else
        at_test = minimum(getfield.(ps, :time))
        if at_test > Time(ereport.par[1])
            _ats = DateTime(Date(ereport.par[1]), at_test) - Hour(8)
        else
            _ats = DateTime(Date(ereport.par[1])+Day(1), at_test) - Hour(8)
        end
        push!(s["phases"], Dict(
            "type" => "S",
            "at" => _ats
        ))
    end
    if isnan(sacf.hdr["scale"])
        s["meta_scale"] = 1.0
    else
        s["meta_scale"] = sacf.hdr["scale"]
    end
    # - Green function setting
    # general options
    s["green_modeltype"] = "3D_COMPRESSED"
    s["green_model"] = "cvm2"
    s["green_tsource"] = hduration(event["magnitude"])
    s["green_dt"] = 0.1

    (glib_station_id, glib_station_dist) = nearest_station_in_glib(s["meta_lat"], s["meta_lon"], glats, glons)
    glib_station_tag = glibindex[glib_station_id][1]
    # baz = SeisTools.Geodesy.azimuth(s["meta_lat"], s["meta_lon"], event["latitude"], event["longitude"])
    n = (event["latitude"] - s["meta_lat"]) * C0
    e = (event["longitude"] - s["meta_lon"]) * C0 * cosd(s["meta_lat"])
    # n = s["base_distance"] * cosd(baz)
    # e = s["base_distance"] * sind(baz)
    s["green_modelpath"] = joinpath(PTH["glib"],
        glib_station_tag,
        Eareacode*".bin"
    )
    # println(s["green_modelpath"])
    (_, _, tp, ts, _) = JuliaSourceMechanism.Green.cglib_readlocation(s["green_modelpath"], n, e, event["depth"])
    # ! DEBUG
    # (_, _, tp, ts, _) = JuliaSourceMechanism.Green.cglib_readlocation(s["green_modelpath"], -170, -150, event["depth"])
    # * phase infomation
    for p in s["phases"]
        # p["tt"] = round(p["at"] - event["origintime"], Millisecond).value * 1e-3
        if p["type"] == "P"
            defaultband = [0.05, 0.2]
            p["tt"] = tp
        else
            defaultband = [0.05, 0.1]
            p["tt"] = ts
        end
        p["xcorr_order"] = 4
        p["xcorr_band"] = defaultband
        p["xcorr_maxlag"] = 1.0 / p["xcorr_band"][2]
        if p["type"] == "P"
            p["xcorr_trim"] = [-2.0, 3.0] ./ p["xcorr_band"][2]
            p["polarity_obs"] = 0.0
            # if isnan(sacf.hdr["a"])
            #     p["polarity_obs"] = 0.0
            # else
            #     shift = round(Int, (sacf.hdr["a"] - sacf.hdr["b"]) / sacf.hdr["delta"]) + 1
            #     p["polarity_obs"] = sign(sum(sacf.data[shift:shift+9]))
            # end
            p["polarity_trim"] = [0.0, s["green_tsource"]]
        else
            p["xcorr_trim"] = [-4.0, 6.0] ./ p["xcorr_band"][2]
            p["polarity_obs"] = NaN
            p["polarity_trim"] = [NaN, NaN]
        end
        local dt = s["meta_dt"]
        tl = p["xcorr_trim"][2] - p["xcorr_trim"][1]
        while (dt + 1e-3) * 200 < tl
            dt += 1e-3
        end
        p["xcorr_dt"] = dt
    end
end

# = autopick arrivals

pickflag = falses(length(stations))
while !all(pickflag)
    i = findfirst(!, pickflag)
    js = findall(s->(s["network"] == stations[i]["network"]) &&
        (s["station"] == stations[i]["station"]), stations)
    sacs = map(js) do j
        open(SeisTools.SAC.read, joinpath(dataroot, "sac", stations[j]["meta_file"]))
    end
    tp = map(sacs) do sacf
        # (it, _) = JuliaSourceMechanism.pick_windowratio(sacf.data, round(Int, 5.0/sacf.hdr["delta"]))
        (it, _) = JuliaSourceMechanism.pick_stalta(sacf.data,
            round(Int, 2.0/sacf.hdr["delta"]), round(Int, 10.0/sacf.hdr["delta"]))
        SeisTools.SAC.DateTime(sacf.hdr) + Millisecond(round(Int, sacf.hdr["b"]*1000 + it*sacf.hdr["delta"]*1000))
    end
    tp = minimum(tp)
    for j in js
        ip = findfirst(p->p["type"] == "P", stations[j]["phases"])
        if abs(stations[j]["phases"][ip]["at"] - TIME_DEFAULT) > Minute(1)
            continue
        end
        shift = round(Int, Millisecond(tp - stations[j]["base_trim"][1])/
            Millisecond(round(Int, stations[j]["meta_dt"]*1000)))
        # l = round(Int, 0.5/stations[j]["meta_dt"])
        # stations[j]["phases"][ip]["polarity_obs"] = sign(sum(sacs[j-minimum(js)+1].data[shift:shift+l]))
        stations[j]["phases"][ip]["polarity_obs"] = 0.0
        stations[j]["phases"][ip]["at"] = tp
    end
    ip = findfirst(p->p["type"] == "P", stations[js[1]]["phases"])
    tp = stations[js[1]]["phases"][ip]["at"]
    w = zeros(minimum(length.(getfield.(sacs, :data))), length(sacs))
    for i = eachindex(sacs)
        w[:, i] .= sacs[i].data[1:size(w, 1)]
    end

    ws = SeisTools.DataProcess.cut(w, SeisTools.SAC.DateTime(sacs[1].hdr), tp-Second(20),
        tp+Second(60), Millisecond(round(Int, sacs[1].hdr["delta"]*1000)))
    wt = deepcopy(ws[2])
    # for i = axes(wt, 2)
    #     wt[:, i] = SeisTools.DataProcess.bandpass(ws[2][:, i], 0.05, 1.0, 1/sacs[1].hdr["delta"])
    # end
    W = round(Int, 10.0/sacs[1].hdr["delta"])
    (_, r1) = JuliaSourceMechanism.pick_windowratio(wt[:, 1], W)
    # (_, r2) = JuliaSourceMechanism.pick_freedom(wt, W)
    (_, it1) = findmax(r1)
    # (_, it2) = findmax(r2)
    it1 += W
    # it2 += W
    for j in js
        is = findfirst(p->p["type"] == "S", stations[j]["phases"])
        if abs(stations[j]["phases"][is]["at"] - TIME_DEFAULT) > Minute(1)
            continue
        end
        stations[j]["phases"][is]["at"] = tp + Millisecond(round(Int, it1*sacs[1].hdr["delta"]*1000))
    end
    for j in js
        pickflag[j] = true
    end
end

# = construct env data

env = Dict("algorithm" => algorithm,
        "event" => event,
        "stations" => stations,
        "dataroot" => dataroot)

JuliaSourceMechanism.loaddata!(env)
for s in env["stations"]
    ot = Millisecond(env["event"]["origintime"] - s["base_begintime"]).value * 1e-3
    (meta, _) = JuliaSourceMechanism.Green.scangreenfile(normpath(env["dataroot"],
                                                                "greenfun",
                                                                @sprintf("%s-%.4f", s["green_model"],
                                                                        env["algorithm"]["searchdepth"]),
                                                                s["network"] * "." * s["station"] * "." *
                                                                s["component"] * ".gf"))
    for p in s["phases"]
        if p["type"] == "P"
            p["tt"] = meta["tp"] + ot
        elseif p["type"] == "S"
            p["tt"] = meta["ts"] + ot
        end
    end
end
status = Dict{String,Any}()
status["saveplotdata"] = true
status["saveplotdatato"] = abspath(".tmpplot.mat")

# = copy data and delete extra station

tenv = Dict{String,Any}()
tenv["algorithm"] = deepcopy(env["algorithm"])
tenv["event"] = deepcopy(env["event"])
tenv["dataroot"] = deepcopy(env["dataroot"])
scpflag = falses(length(env["stations"]))
cpcheck = falses(length(env["stations"]))
while !all(scpflag)
    i = findfirst(!, scpflag)
    js = findall(s->(s["network"] == env["stations"][i]["network"]) &&
        (s["station"] == env["stations"][i]["station"]), env["stations"])
    w = 0.0
    for j in js
        for p in env["stations"][j]["phases"]
            w += ("xcorr_weight" in keys(p)) ? p["xcorr_weight"] : 1.0
        end
    end
    if !iszero(w)
        for j in js
            cpcheck[j] = true
        end
    end
    for j in js
        scpflag[j] = true
    end
end
tenv["stations"] = Vector{Dict{String,Any}}(undef, count(cpcheck))
cpid = 1
for i = eachindex(env["stations"])
    global cpid
    if cpcheck[i]
        tenv["stations"][cpid] = deepcopy(env["stations"][i])
        cpid += 1
    end
end

jldsave(joinpath(dataroot, "auto.jld2"); env=tenv, status=status)
