# using JuliaSourceMechanism, JLD2, LinearAlgebra, Printf, Dates, SeisTools, Statistics

Geo = SeisTools.Geodesy
SS = SeisTools.Source
SDP = SeisTools.DataProcess

msecond(x::Real) = Millisecond(round(Int, x * 1000.0))

function qualitycontrol(env)
    nenv = Setting()
    nenv["dataroot"] = env["dataroot"]
    nenv["algorithm"] = env["algorithm"]
    nenv["event"] = env["event"]
    nenv["stations"] = Setting[]
    for s in env["stations"]
        flag = true
        flag &= !SeisTools.QualityControl.islimitedamplitude(s["base_record"], round(Int, 1.0 / s["meta_dt"]))
        # flag &= !SeisTools.QualityControl.hasconstrecord(s["base_record"], round(Int, 5.0 / s["meta_dt"]))
        w = 0.0
        for p in s["phases"]
            if haskey(p, "xcorr_weight")
                w += p["xcorr_weight"]
            else
                w += 1.0
            end
        end
        flag &= w > 0.0
        if flag
            push!(nenv["stations"], deepcopy(s))
        end
    end
    return nenv
end

function loadtp!(env)
    for s in env["stations"]
        ot = Millisecond(env["event"]["origintime"] - s["base_begintime"]).value * 1e-3
        (meta, _) = Green.scangreenfile(normpath(env["dataroot"], "greenfun",
            @sprintf("%s-%.4f", s["green_model"], env["algorithm"]["searchdepth"]),
            s["network"] * "." * s["station"] * "." * s["component"] * ".gf"))
        for p in s["phases"]
            if p["type"] == "P"
                p["tt"] = meta["tp"] + ot
            elseif p["type"] == "S"
                p["tt"] = meta["ts"] + ot
            end
        end
    end
end

timediff_sec(x::DateTime, y::DateTime) = Millisecond(x - y).value * 1e-3

function inverse_focalmech!(nenv, misfits)
    loadtp!(nenv)
    preprocess!(nenv, misfits)
    _STEPS = [(11, 13, 14), (5, 3, 3), (1, 1, 1)]

    Grid.setstart!(0.0, 0.0, -90.0)
    Grid.setstep!(_STEPS[1][1], _STEPS[1][2], _STEPS[1][3])
    if iszero(mod(360, _STEPS[1][1]))
        Grid.setstop!(360.0 - _STEPS[1][1], 90.0, 90.0)
    else
        Grid.setstop!(floor(360.0/_STEPS[1][1])*_STEPS[1][1], 90.0, 90.0)
    end

    mech = zeros(3)
    mech0 = zeros(3, length(_STEPS))
    minval = 0.0
    minval_xcorr = 0.0
    minval_pol = 0.0

    for it in eachindex(_STEPS)
        # global mech, minval, minval_xcorr, minval_pol
        (sdr, phaselist, _, misfitdetail) = inverse!(nenv, misfits, Grid)
        weight = normalize(map(x -> x[1].weight(x[2], nenv, nenv), phaselist), 1)
        weight_xcorr = normalize(map(x->x[1] == XCorr ? 1.0 : 0.0, phaselist), 1)
        weight_pol = normalize(map(x->x[1] == Polarity ? 1.0 : 0.0, phaselist), 1)
        totalmisfit = replace(misfitdetail, NaN => 0.0) * weight
        mis_xcorr = replace(misfitdetail, NaN => 0.0) * weight_xcorr
        mis_pol = replace(misfitdetail, NaN => 0.0) * weight_pol
        (minval, minidx) = findmin(totalmisfit)
        minval_xcorr = mis_xcorr[minidx]
        minval_pol = mis_pol[minidx]
        mech .= sdr[minidx]
        mech0[:, it] .= mech
        if it < length(_STEPS)
            Grid.setstart!(mech[1]-2*_STEPS[it][1],
                max(0.0, mech[2]-2*_STEPS[it][2]),
                max(-90.0, mech[3]-2*_STEPS[it][3]))
            Grid.setstep!(_STEPS[it+1]...)
            Grid.setstop!(mech[1]+2*_STEPS[it][1],
                min(90.0, mech[2]+2*_STEPS[it][2]),
                min(90.0, mech[3]+2*_STEPS[it][3]))
        end
    end
    for i = 1:length(_STEPS)-1
        flag = true
        _ag = mod(abs(mech0[1, i+1] - mech0[1, i]), 360.0)
        flag &= min(_ag, 360.0-_ag) <= 2 * _STEPS[i][1]
        flag &= abs(mech0[2, i+1] - mech0[2, i]) <= 2 * _STEPS[i][2]
        flag &= abs(mech0[3, i+1] - mech0[3, i]) <= 2 * _STEPS[i][3]
        if !flag
            printstyled(mech0[:, i], " -> ", mech0[:, i+1], "step:", _STEPS[i], "\n", color=:yellow)
        end
    end
    return (mech, minval, minval_xcorr, minval_pol)
end

function inverse_momenttensor!(mt_final, nenv, thres::Real, nit::Integer=2)
    for _ in 1:nit
        maxcorrval = Float64[]
        rec_seg = Vector{Float64}[]
        g_seg = Matrix{Float64}[]

        for ista in eachindex(nenv["stations"])
            s = nenv["stations"][ista]
            for p in s["phases"]
                pdt = Millisecond(round(Int, p["xcorr_dt"] * 1e3))
                w_resample = SDP.resample(s["base_record"], s["meta_dt"], p["xcorr_dt"])
                w_filt = SDP.bandpass(w_resample, p["xcorr_band"][1], p["xcorr_band"][2], 1.0/p["xcorr_dt"])
                w_cut = SDP.cut(w_filt, s["base_begintime"],
                    p["at"] + Millisecond(round(Int, p["xcorr_trim"][1] * 1e3)),
                    p["at"] + Millisecond(round(Int, p["xcorr_trim"][2] * 1e3)), pdt)
                _rec = w_cut[2]

                g_resample = SDP.resample(s["green_fun"], s["green_dt"], p["xcorr_dt"])
                g_filt = SDP.bandpass(g_resample, p["xcorr_band"][1], p["xcorr_band"][2], 1.0/p["xcorr_dt"])
                g_cut = SDP.cut(g_filt, round(Int, (p["tt"]+p["xcorr_trim"][1])/p["xcorr_dt"]), length(_rec))
                _syn = g_cut * mt_final
                mlag = round(Int, p["xcorr_maxlag"]/p["xcorr_dt"])
                xc = SDP.xcorr_t(_rec, _syn, -mlag, mlag)
                (vmax, imax) = findmax(vec(xc.c))
                push!(maxcorrval, vmax/norm(_rec)/norm(_syn))
                shift = xc.lag[imax]
                L = length(_rec)
                push!(g_seg, p["xcorr_greenfun"][1+max(0, -shift):L+min(0, -shift), :])
                push!(rec_seg, _rec[1+max(0, shift):L+min(0, shift)])
            end
        end

        sp = sortperm(maxcorrval, rev=true)
        ncut = findlast(>=(thres), maxcorrval[sp])
        ncut = isnothing(ncut) ? length(sp) : ncut
        ncut = (ncut < 6) ? length(sp) : ncut

        G = zeros(0, 6)
        d = zeros(0)

        for i = 1:ncut
            # A = maximum(abs, rec_seg[sp[i]])
            A = sqrt(sum(abs2, g_seg[sp[i]])/size(g_seg[sp[i]], 1))
            G = vcat(G, g_seg[sp[i]]./A)
            append!(d, rec_seg[sp[i]]./A)
        end

        total_mt = (G'*G)\(G'*d)
        n1 = sum(abs2, total_mt) + sum(abs2, total_mt[4:6])
        total_mt ./= sqrt(n1)/sqrt(2)
        dmt = total_mt - mt_final
        mt_final .= total_mt
        if (sum(abs2, dmt) + sum(abs2, dmt[4:6])) < 0.0001 * 2.0
            break
        end
    end
end

function inverse_tsource!(tslist::Vector{<:Real}, nenv)
    mechlist = zeros(3, length(tslist))
    minvallist = zeros(length(tslist))
    minval_xcorrlist = zeros(length(tslist))
    minval_pollist = zeros(length(tslist))

    for i = eachindex(tslist)
        for s in nenv["stations"]
            s["green_tsource"] = tslist[i]
        end
        preprocess!(nenv, misfits)
        (mechlist[:, i], minvallist[i], minval_xcorrlist[i], minval_pollist[i]) = inverse_focalmech!(nenv, misfits)
    end

    imin = argmin(minvallist)
    return tslist[imin]
end

function inverse_depth(deplist::AbstractVector{<:Real}, nenv, misfits)
    val = zeros(length(deplist))
    for i = eachindex(deplist)
        # println("depth: ", deplist[i])
        tenv = deepcopy(nenv)
        tenv["algorithm"]["searchdepth"] = deplist[i]
        JuliaSourceMechanism.calcgreen!(tenv)
        (_, val[i], _, _) = inverse_focalmech!(tenv, misfits)
    end
    return val
end

function reselect_channel(env, mt, thres::Real)
    nenv = Setting()
    nenv["dataroot"] = env["dataroot"]
    nenv["algorithm"] = env["algorithm"]
    nenv["event"] = env["event"]
    nenv["stations"] = Setting[]
    for s in env["stations"]
        t = Setting()
        for k in keys(s)
            if k == "phases"
                continue
            end
            t[k] = deepcopy(s[k])
        end
        ip = findfirst(_p->_p["type"]=="P", s["phases"])
        if !isnothing(ip)
            p = s["phases"][ip]
            _rec = p["xcorr_record"]
            _syn = p["xcorr_greenfun"] * mt
            mlag = round(Int, p["xcorr_maxlag"] / p["xcorr_dt"])
            xc = SDP.xcorr_t(_rec, _syn, -mlag, mlag)
            maxc1 = maximum(xc.c)/norm(_rec)/norm(_syn)
            if maxc1 >= thres
                if !haskey(t, "phases")
                    t["phases"] = Setting[]
                end
                push!(t["phases"], deepcopy(s["phases"][ip]))
            end
        else
            maxc1 = 0.0
        end

        is = findfirst(_p->_p["type"]=="S", s["phases"])
        if !isnothing(is)
            p = s["phases"][is]
            _rec = p["xcorr_record"]
            _syn = p["xcorr_greenfun"] * mt
            mlag = round(Int, p["xcorr_maxlag"] / p["xcorr_dt"])
            xc = SDP.xcorr_t(_rec, _syn, -mlag, mlag)
            maxc2 = maximum(xc.c)/norm(_rec)/norm(_syn)
            if maxc2 >= thres
                if !haskey(t, "phases")
                    t["phases"] = Setting[]
                end
                push!(t["phases"], deepcopy(s["phases"][is]))
            end
        else
            maxc2 = 0.0
        end

        if (maxc1 >= thres) || (maxc2 >= thres)
            push!(nenv["stations"], t)
        end
    end
    return nenv
end

function correctpick!(env, refmt, pname)
    needset = trues(length(env["stations"]))
    while any(needset)
        i = findfirst(needset)
        js = findall(env["stations"]) do _s
                (_s["network"] == env["stations"][i]["network"]) &&
                (_s["station"] == env["stations"][i]["station"])
            end
        juse = findall(env["stations"]) do _s
                (_s["network"] == env["stations"][i]["network"]) &&
                (_s["station"] == env["stations"][i]["station"]) &&
                (pname in getindex.(_s["phases"], "type"))
            end
        if isempty(juse)
            for j in js
                needset[j] = false
            end
            continue
        end
        gdt = env["stations"][i]["green_dt"]
        winN = size(env["stations"][i]["green_fun"], 1)
        w = zeros(winN, length(juse))
        g = zeros(winN, length(juse))
        for j = eachindex(juse)
            SDP.resample!(@view(w[:, j]), env["stations"][juse[j]]["base_record"] |> SDP.detrend |> SDP.taper)
            g[:, j] .=  env["stations"][juse[j]]["green_fun"] * refmt
        end
        iip = let
            ips_ = map(env["stations"][juse]) do _s
                findfirst(_p->_p["type"]==pname, _s["phases"])
            end |> unique
            all(isnothing, ips_) ? nothing : ips_[findfirst(!isnothing, ips_)]
        end
        if isnothing(iip)
            for j in js
                needset[j] = false
            end
            continue
        end
        pAT = env["stations"][i]["phases"][iip]["at"]
        pTT = env["stations"][i]["phases"][iip]["tt"]
        (f1, f2) = env["stations"][i]["phases"][iip]["xcorr_band"]
        (lwin, rwin) = env["stations"][i]["phases"][iip]["xcorr_trim"]
        pWinN = round(Int, (rwin - lwin) / f0 / gdt)
        prshift = round(Int, (timediff_sec(pAT, env["stations"][i]["base_begintime"]) + lwin / f0) / gdt)
        pShiftN = round(Int, (pTT + lwin / f0)/gdt)
        wt = SDP.bandpass(w, f1, f2, 1.0/gdt)
        gt = SDP.bandpass(g, f1, f2, 1.0/gdt)
        _w = SDP.cut(wt, prshift, pWinN)
        _g = SDP.cut(gt, pShiftN, pWinN)
        xcvalue = zeros(size(_w, 1)+size(_g, 1)-1, length(juse))
        for j = eachindex(juse)
            xc = SDP.xcorr_t(_w[:, j], _g[:, j])
            xcvalue[:, j] .= xc.c ./ norm(_w[:, j]) ./ norm(_g[:, j])
        end
        xcv = sum(xcvalue, dims=2) |> x->reshape(x, length(x))
        (_, maxi) = findmax(xcv)
        psft = (maxi - size(_g, 1))*gdt
        # psft = pshift[imaxfreqp]

        for j in juse
            ip = findfirst(_p->_p["type"] == pname, env["stations"][j]["phases"])
            if isnothing(ip)
                continue
            end
            env["stations"][j]["phases"][ip]["at"] -= Millisecond(round(Int, psft*1000))
        end

        for j in js
            needset[j] = false
        end
    end
    return env
end

function update_filter_band!(env, refmt::Vector{<:Real}, ftrylist::AbstractVector{<:Real}, pname::String,
        lwin::Real, rwin::Real)
    needset = trues(length(env["stations"]))
    while any(needset)
        i = findfirst(needset)
        js = findall(env["stations"]) do _s
                (_s["network"] == env["stations"][i]["network"]) &&
                (_s["station"] == env["stations"][i]["station"])
            end
        juse = findall(env["stations"]) do _s
                (_s["network"] == env["stations"][i]["network"]) &&
                (_s["station"] == env["stations"][i]["station"]) &&
                (pname in getindex.(_s["phases"], "type"))
            end
        if isempty(juse)
            for j in js
                needset[j] = false
            end
            continue
        end
        gdt = env["stations"][i]["green_dt"]
        winN = size(env["stations"][i]["green_fun"], 1)
        w = zeros(winN, length(juse))
        g = zeros(winN, length(juse))
        # println("refmt:", refmt, ", winN: ", winN)
        for j = eachindex(juse)
            SDP.resample!(@view(w[:, j]), env["stations"][juse[j]]["base_record"] |> SDP.detrend |> SDP.taper)
            # println(size(g), " ", size(env["stations"][juse[j]]["green_fun"]))
            g[:, j] .=  env["stations"][juse[j]]["green_fun"] * refmt
        end
        # iip = findfirst(_p->_p["type"]==pname, env["stations"][i]["phases"])
        iip = let
            ips_ = map(env["stations"][juse]) do _s
                findfirst(_p->_p["type"]==pname, _s["phases"])
            end |> unique
            all(isnothing, ips_) ? nothing : ips_[findfirst(!isnothing, ips_)]
        end
        if isnothing(iip)
            for j in js
                needset[j] = false
            end
            continue
        end
        pAT = env["stations"][juse[1]]["phases"][iip]["at"]
        pTT = env["stations"][juse[1]]["phases"][iip]["tt"]
        pmax_xc = zeros(length(ftrylist))
        pshift = zeros(length(ftrylist))
        Threads.@threads for ifreq = eachindex(ftrylist)
            local f0 = ftrylist[ifreq]
            pWinN = round(Int, (rwin - lwin) / f0 / gdt)
            prshift = round(Int, (timediff_sec(pAT, env["stations"][i]["base_begintime"]) + lwin / f0) / gdt)
            pShiftN = round(Int, (pTT + lwin / f0)/gdt)
            wt = SDP.bandpass(w, 0.05, f0, 1.0/gdt)
            gt = SDP.bandpass(g, 0.05, f0, 1.0/gdt)
            local _w = SDP.cut(wt, prshift, pWinN)
            local _g = SDP.cut(gt, pShiftN, pWinN)
            local xcvalue = zeros(size(_w, 1)+size(_g, 1)-1, length(juse))
            for j = eachindex(juse)
                xc = SDP.xcorr_t(_w[:, j], _g[:, j])
                xcvalue[:, j] .= xc.c ./ norm(_w[:, j]) ./ norm(_g[:, j])
            end
            xcv = sum(xcvalue, dims=2) |> x->reshape(x, length(x))
            (maxv, maxi) = findmax(xcv)
            pmax_xc[ifreq] = maxv
            pshift[ifreq] = (maxi - size(_g, 1))*gdt
        end
        (_, imaxfreqp) = findmax(pmax_xc)
        pf0 = ftrylist[imaxfreqp]
        # psft = pshift[imaxfreqp]

        for j in juse
            ip = findfirst(_p->_p["type"] == pname, env["stations"][j]["phases"])
            if isnothing(ip)
                continue
            end
            # env["stations"][j]["phases"][ip]["at"] += Millisecond(round(Int, psft*1000))
            env["stations"][j]["phases"][ip]["xcorr_band"] = [0.05, pf0]
            env["stations"][j]["phases"][ip]["xcorr_trim"] = [lwin, rwin] ./ pf0
            env["stations"][j]["phases"][ip]["xcorr_maxlag"] = 1.0 / pf0
            env["stations"][j]["phases"][ip]["xcorr_dt"] = 1 / pf0 / 20
        end

        for j in js
            needset[j] = false
        end
    end
    return env
end

function locatepoint(x::Real, xs::AbstractVector{<:Real}, ERRoutofrange::Bool=true)
    if ERRoutofrange && ((x < minimum(xs)) || (x > maximum(xs)))
        error(string(x)*" Out of Model Range")
    end
    i = findfirst(>(x), xs)
    if isnothing(i)
        i = length(xs)
        h = 1.0
    elseif i == 1
        i = 2
        h = 0.0
    else
        h = (x - xs[i-1]) / (xs[i] - xs[i-1])
    end
    return (i, h)
end

# function bilinear(x, i, j, p, q)
#     return x[i-1, j-1] * (1.0 - p) * (1.0 - q) +
#            x[i-1, j] * (1.0 - p) * q +
#            x[i, j-1] * p * (1.0 - q) +
#            x[i, j] * p * q
# end

function triplelinear(x, i, j, k, p, q, r)
    if any(iszero, x[i-1:i, j-1:j, k-1:k])
        error("Velocity is Zero")
    end
    return x[i-1, j-1, k-1] * (1.0 - p) * (1.0 - q) * (1.0 - r) +
           x[i-1, j-1, k] * (1.0 - p) * (1.0 - q) * r +
           x[i-1, j, k-1] * (1.0 - p) * q * (1.0 - r) +
           x[i-1, j, k] * (1.0 - p) * q * r +
           x[i, j-1, k-1] * p * (1.0 - q) * (1.0 - r) +
           x[i, j-1, k] * p * (1.0 - q) * r +
           x[i, j, k-1] * p * q * (1.0 - r) +
           x[i, j, k] * p * q * r
end

function calmag(env, mt, mag_a, mag_b)
    (latlist, lonlist, deplist, mvp, mvs) = let
        io = open(joinpath(@__DIR__, "SWChinaCVM2.0.sea_level.bin"))
        nlat = read(io, Int32)
        nlon = read(io, Int32)
        ndep = read(io, Int32)
        lats = zeros(Float32, nlat)
        lons = zeros(Float32, nlon)
        deps = zeros(Float32, ndep)
        read!(io, lats)
        read!(io, lons)
        read!(io, deps)
        vp = zeros(Float32, ndep, nlon, nlat)
        vs = zeros(Float32, ndep, nlon, nlat)
        read!(io, vp)
        read!(io, vs)
        close(io)
        (lats, lons, deps, vp, vs)
    end

    (mlat, p) = locatepoint(env["event"]["latitude"], latlist)
    (mlon, q) = locatepoint(env["event"]["longitude"], lonlist)
    (mdep, r) = locatepoint(env["algorithm"]["searchdepth"], deplist, false)
    srcvs = triplelinear(mvs, mdep, mlon, mlat, r, q, p)
    E0 = 2.7*srcvs*srcvs*1e15

    wrs = Vector{Float64}[]
    wss = Vector{Float64}[]
    gAs = Float64[]
    A0 = Float64[]
    for is = eachindex(env["stations"])
        s = env["stations"][is]
        # stag = s["network"]*"."*s["station"]
        # sac = SeisTools.SAC.read(joinpath("../../sac_rr/", s["meta_file"]))
        # _w = SDP.cut(sac.data, SeisTools.SAC.DateTime(sac.hdr),
        #     s["base_trim"][1], s["base_trim"][2], _Second(s["meta_dt"]))
        _w = s["base_record"]
        SDP.detrend!(_w)
        SDP.taper!(_w)
        _w_resample = zeros(size(s["green_fun"], 1))
        SDP.resample!(_w_resample, _w)
        gA = sqrt(sum(abs2, s["green_fun"])/size(s["green_fun"], 1))
        for p = s["phases"]
            _wf = SDP.bandpass(_w_resample,
                p["xcorr_band"][1], p["xcorr_band"][2], 1/s["green_dt"])
            at_r = round(Int, (Millisecond(p["at"] - s["base_trim"][1])+
                msecond(p["xcorr_trim"][1]))/msecond(s["green_dt"]))
            Nmaxlag = round(Int, p["xcorr_maxlag"] / p["xcorr_dt"])
            _corr = SDP.xcorr_t(p["xcorr_record"], p["xcorr_greenfun"]*mt, -Nmaxlag, Nmaxlag)
            sshift = _corr.lag[argmax(vec(_corr.c))] * p["xcorr_dt"]
            at_s = round(Int, (p["tt"]+sshift+p["xcorr_trim"][1])/s["green_dt"])
            ntrim = round(Int, (p["xcorr_trim"][2]-p["xcorr_trim"][1])/s["green_dt"])
            wr_seg = zeros(ntrim)
            ws_seg = zeros(ntrim)
            SDP.cut!(wr_seg, _wf, at_r, 1, ntrim)
            _g = SDP.bandpass(s["green_fun"]*mt,
                p["xcorr_band"][1], p["xcorr_band"][2], 1/s["green_dt"])
            SDP.cut!(ws_seg, _g, at_s, 1, ntrim)
            corr = SDP.xcorr_t(wr_seg, ws_seg)
            sft2 = corr.lag[argmax(corr.c)]
            SDP.cut!(ws_seg, _g, at_s-sft2, 1, ntrim)
            push!(wrs, wr_seg)
            push!(wss, ws_seg)
            push!(gAs, gA)
            push!(A0, s["meta_scale"])
        end
    end

    Mtest = zeros(length(wss))
    for i = eachindex(Mtest)
        rmIdxs = [1:i-1; i+1:length(wss)]
        _wr = Float64[]
        _ws = Float64[]
        for j = eachindex(wrs)
            if j in rmIdxs
                continue
            end
            append!(_wr, wrs[j] ./ (A0[j] * gAs[j]))
            append!(_ws, wss[j] ./ gAs[j])
        end
        Mtest[i] = sum(_wr .* _ws) / sum(abs2, _ws)
    end
    M0 = E0*mean(Mtest)*1e10
    # return (mean(Mtest), mag_a*log10(mean(Mtest)) + mag_b)
    # return (M0, (2//3)*log10(M0)-15.75)
    return (M0, 0.6027566*log10(M0)-20.9887342)
end


function channeltag(s::String)
    _l = split(s, '.')
    return String(_l[2] * _l[1] * _l[3])
end

function resptag(s::String)
    _l = split(s, '.')
    return String("RESP."*join(_l[1:4], '.'))
end

const charlist = let
    t = String[]
    for i = 1:26, j = 1:26
        push!(t, String([Char(i+64), Char(j+64)]))
    end
    Tuple(t)
end

lat_tag(x::Real) = charlist[round(Int, 2*(x+0.25) + 180.0)]

function nearest_center(x::Real)
    a = floor((x-0.25)*2)*0.5+0.25
    b = ceil((x-0.25)*2)*0.5+0.25
    if abs(x-a) < abs(x-b)
        return a
    else
        return b
    end
end

function area_code(lat::Real, lon::Real)
    return charlist[round(Int, 2*(nearest_center(lat)+0.25) + 180.0)]*
        @sprintf("%03d", round(Int, (nearest_center(lon)+0.25)*2.0))
end

function nearest_station_in_glib(lat::Real, lon::Real, glat, glon)
    dvec = map(eachindex(glat)) do i
        Geo.distance(lat, lon, glat[i], glon[i])
    end
    mi = argmin(dvec)
    return (mi, dvec[mi]*1e-3)
end

function hduration(m::Real)
    return exp(0.75734675*m-3.87327442)
end


macro write_result(fname)
    return quote
        mt_final = dc2ts(mech)
        inverse_momenttensor!(mt_final, nenv, 0.7, 1)

        result = Setting()
        for s in nenv["stations"]
            sta = s["network"] * "." * s["station"]
            if !(sta ∈ keys(result))
                result[sta] = Setting()
            end
            c = s["component"]
            result[sta]["dist"] = s["base_distance"]
            result[sta][c] = Setting()
            for p in s["phases"]
                t = Setting()
                if "xcorr_weight" in keys(p)
                    t["xcorr_weight"] = p["xcorr_weight"]
                else
                    t["xcorr_weight"] = 1.0
                end
                t["xcorr_rec"] = normalize(p["xcorr_record"], Inf)
                t["xcorr_syn"] = normalize(p["xcorr_greenfun"] * dc2ts(mech), Inf)
                t["mt_syn"] = normalize(p["xcorr_greenfun"] * mt_final, Inf)
                mlag = round(Int, p["xcorr_maxlag"] / p["xcorr_dt"])
                xc = SDP.xcorr_t(t["xcorr_rec"], t["mt_syn"], -mlag, mlag)
                t["mt_shift"] = xc.lag[argmax(vec(xc.c))] * p["xcorr_dt"]
                t["xcorr_shift"] = XCorr.detail(p, dc2ts(mech))
                t["xcorr_dt"] = p["xcorr_dt"]
                t["polarity_rec"] = p["polarity_obs"]
                if p["type"] == "P"
                    t["polarity_syn"] = sign(sum(p["polarity_syn"] .* dc2ts(mech)))
                end
                result[sta][c][p["type"]] = t
            end
        end

        (M0, mag) = calmag(nenv, mt_final, 1.0, 0.0)
        result["info_mech"] = mech
        result["info_mt"] = mt_final
        result["info_M0diff"] = M0
        result["info_mag"] = mag
        result["info_misfit"] = minval
        result["info_misvsdep"] = (hlist, misfitlist)
        result["info_misfit_xcorr"] = minval_xcorr
        result["info_misfit_pol"] = minval_pol

        _MT = SeisTools.Source.MomentTensor(mech...)
        pln = SeisTools.Source.focalmechanism(_MT)
        if !isdir(joinpath(dataroot, "result"))
            mkdir(joinpath(dataroot, "result"))
        end
        jldsave(joinpath(dataroot, "result", $(esc(fname*".jld2"))); env=nenv, status=status, result=result)
        open(joinpath(dataroot, "result", $(esc(fname*".txt"))), "w") do io
            @printf(io, "%04d-%02d-%02d %02d:%02d:%02d %.4f %.4f %.4f %.1f %d %d %d %d %d %d\n",
                year(nenv["event"]["origintime"]), month(nenv["event"]["origintime"]),
                day(nenv["event"]["origintime"]), hour(nenv["event"]["origintime"]),
                minute(nenv["event"]["origintime"]), second(nenv["event"]["origintime"]),
                nenv["event"]["latitude"], nenv["event"]["longitude"],
                nenv["algorithm"]["searchdepth"], mag,
                round(Int, pln.plane1[1]), round(Int, pln.plane1[2]), round(Int, pln.plane1[3]),
                round(Int, pln.plane2[1]), round(Int, pln.plane2[2]), round(Int, pln.plane2[3]))
        end
    end
end
