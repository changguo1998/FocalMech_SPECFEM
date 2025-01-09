using DelimitedFiles, SeisTools, TOML, Dates, JuliaSourceMechanism,
    SourceMechUI, Printf, JLD2, Statistics, LinearAlgebra

Geo = SeisTools.Geodesy

include("../env.jl")

include(joinpath(PTH["invsrc"], "lib/io.jl"))
include(joinpath(PTH["invsrc"], "lib/lib.jl"))
include(joinpath(PTH["invscript"], "multistage_lib.jl"))
include("/home/guochang/Projects/NLLOCshell/NLLOCshell.jl")
include("/home/guochang/Projects/glibSEMshell/include/io.jl")

eventname = ARGS[1]

(env, status) = let
    t = load(joinpath(PTH["dat_selectstation"], eventname, "auto.jld2"))
    (t["env"], t["status"])
end

# @info "Copy data"
rawenv = Setting()
rawenv["dataroot"] = abspath(PTH["dat_selectstation"], eventname)
dataroot=rawenv["dataroot"]
rawenv["algorithm"] = env["algorithm"]
rawenv["event"] = env["event"]
rawenv["stations"] = Setting[]
for s in env["stations"]
    t = deepcopy(s)
    push!(rawenv["stations"], t)
end
# @info "Change parameter"
# rawenv["algorithm"]["misfit"] = ["xcorr", "pol"]
rawenv["algorithm"]["weight"] = [1.0, 0.4]

# rawenv = qualitycontrol(rawenv)

misfits = Module[]
for m in rawenv["algorithm"]["misfit"], f in [XCorr, Polarity, PSR, DTW, AbsShift, RelShift]
    if m in f.tags
        push!(misfits, f)
    end
end

# @info "Run"
JuliaSourceMechanism.calcgreen!(rawenv)

(mech, minval, minval_xcorr, minval_pol) = inverse_focalmech!(rawenv, misfits)

nenv = deepcopy(rawenv)
hlist = [nenv["algorithm"]["searchdepth"]]
misfitlist = [minval]
@write_result "result_stage1"

(mech, minval, minval_xcorr, _minval_pol) = let
    local _tfms = zeros(3, 0)
    local _mech = zeros(3)
    local _tmech = zeros(3)
    _mech .= mech
    local _minval = 0.0
    local _minval_xcorr = 0.0
    local _minval_pol = 0.0
    local Nsample = 3
    local Ncount = 0
    while true
        _rmech = randn(3)
        _rmech[1] = min(3.0, max(-3.0, _rmech[1]))
        _rmech[2] = min(3.0, max(-3.0, _rmech[2]))
        _rmech[3] = min(3.0, max(-3.0, _rmech[3]))
        _tmech[1] = mod(_mech[1]+_rmech[1], 360.0)
        _tmech[2] = max(0.0, min(90.0, _mech[2]+_rmech[2]))
        _tmech[3] = min(90.0, max(-90.0, _mech[3] + _rmech[3]))
        update_filter_band!(nenv, dc2ts(_tmech), 0.2:0.01:0.25, "P", -2.0, 3.0)
        update_filter_band!(nenv, dc2ts(_tmech), 0.1:0.01:0.15, "S", -4.0, 6.0)
        (_mech, _minval, _minval_xcorr, _minval_pol) = inverse_focalmech!(nenv, misfits)
        _tfms = hcat(_tfms, _mech)
        if size(_tfms, 2) < Nsample
            # println(mech)
            continue
        end
        cstrike = cosd.(_tfms[1,end-Nsample+1:end])
        sstrike = sind.(_tfms[1,end-Nsample+1:end])
        std1 = sqrt(std(cstrike)^2+std(sstrike)^2)
        std2 = std(_tfms[2,end-Nsample+1:end])
        std3 = std(_tfms[3,end-Nsample+1:end])
        # println([mech; std1; std2; std3])
        if std1 < sind(5) && max(std2, std3) < 2.5
            break
        end
        Ncount += 1
        if Ncount > 5
            break
        end
    end
    (_mech, _minval, _minval_xcorr, _minval_pol)
end

h = nenv["algorithm"]["searchdepth"]
r = 5.0
dh = 1.0
hlist = Float64[]
misfitlist = Float64[]
for _ = 1:10
    global h
    local tl = filter(_h->!(_h in hlist), max(0.0, h-r):dh:min(50.0, h+r))
    append!(hlist, tl)
    local val = inverse_depth(tl, nenv, misfits)
    append!(misfitlist, val)
    local h2 = hlist[argmin(misfitlist)]
    if abs(h2 - h) < 2.0
        h = h2
        break
    else
        h = h2
    end
end

nenv["algorithm"]["searchdepth"] = h
(mech, minval, minval_xcorr, minval_pol) = inverse_focalmech!(nenv, misfits)

@write_result "result_stage2"

nenv = reselect_channel(nenv, dc2ts(mech), 0.6)
(mech, minval, minval_xcorr, minval_pol) = inverse_focalmech!(nenv, misfits)

h = nenv["algorithm"]["searchdepth"]
r = 5.0
dh = 1.0
hlist = Float64[]
misfitlist = Float64[]
for _ = 1:10
    global h
    local tl = filter(_h->!(_h in hlist), max(0.0, h-r):dh:min(50.0, h+r))
    append!(hlist, tl)
    local val = inverse_depth(tl, nenv, misfits)
    append!(misfitlist, val)
    local h2 = hlist[argmin(misfitlist)]
    if abs(h2 - h) < 2.0
        h = h2
        break
    else
        h = h2
    end
end

nenv["algorithm"]["searchdepth"] = h
(mech, minval, minval_xcorr, minval_pol) = inverse_focalmech!(nenv, misfits)

@write_result "result_stage3"
