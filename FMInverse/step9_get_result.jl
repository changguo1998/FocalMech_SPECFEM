using DelimitedFiles, TOML, JLD2, SeisTools, Printf, Dates, LinearAlgebra

include("env.jl")
include(joinpath(PTH["invsrc"], "lib/io.jl"))
include(joinpath(PTH["invsrc"], "lib/lib.jl"))

function datestr(t::DateTime)
    return @sprintf("%04d-%02d-%02d %02d:%02d:%02d.%03d", (t .|> [year, month, day, hour, minute, second, millisecond])...)
end

function maxdelta(azs)
    t1 = sort(azs)
    t2 = diff(t1)
    push!(t2, t1[1]-t1[end]+360)
    return maximum(t2)
end

function networkratio(x, y)
    S = svd(permutedims([x y]))
    return S.S[2]/S.S[1]
end

function getinfo(e::AbstractString)
    edir = joinpath(datdir, e)
    rdat = JLD2.load(joinpath(edir, "result", "result_stage3.jld2"))
    evt = rdat["env"]["event"]
    res = rdat["result"]
    stas = rdat["env"]["stations"]
    geos = Dict()
    for s in stas
        tag = join([s["network"], s["station"]], '.')
        if tag ∉ keys(geos)
            geos[tag] = Dict("dist"=>s["base_distance"],
                "az"=>s["base_azimuth"])
        end
    end

    strbuf = Vector{String}(undef, 13)
    strbuf[1] = datestr(evt["origintime"])
    strbuf[2] = @sprintf("%.3f", evt["latitude"])
    strbuf[3] = @sprintf("%.3f", evt["longitude"])
    (d, m) = res["info_misvsdep"]
    strbuf[4] = @sprintf("%.1f", d[argmin(m)])
    strbuf[5] = @sprintf("%.3e", res["info_M0diff"])
    strbuf[6] = @sprintf("%.1f", res["info_mag"])
    strbuf[7] = @sprintf("%3d", res["info_mech"][1])
    strbuf[8] = @sprintf("%2d", res["info_mech"][2])
    strbuf[9] = @sprintf("%4d", res["info_mech"][3])
    strbuf[10] = @sprintf("%.6f", minimum(m))

    rsta = filter(!startswith("info_"), collect(keys(res)))

    strbuf[11] = @sprintf("%3d", length(rsta))

    dists = map(k->geos[k]["dist"], rsta)
    azs = map(k->geos[k]["az"], rsta)

    strbuf[12] = @sprintf("%3d", round(Int, maxdelta(azs)))

    x = dists .* cosd.(azs)
    y = dists .* sind.(azs)

    strbuf[13] = length(rsta) > 3 ? @sprintf("%.6f", networkratio(x, y)) : "0.0"
    return join(strbuf, ' ')
end

datdir = PTH["dat_selectstation"]

events = filter(d->isfile(joinpath(datdir, d, "result/result_stage3.jld2")), readdir(datdir))

open(io->foreach(l->println(io, l), map(getinfo, events)), "result_e.txt", "w")
