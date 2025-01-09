using TOML, Dates

function locategreenlibid(stationpath::AbstractString, dist::Real, az::Real)
    s = TOML.parsefile(joinpath(stationpath, "setting.toml"))
    x = dist * cosd(az)
    y = dist * sind(az)
    n = x + s["srcloc"][1]
    e = y + s["srcloc"][2]
    ids = Int[]
    for r = s["receiver"]
        if (r["n"][1] <= n) && (n <= r["n"][3]) && (r["e"][1] <= e) && (e <= r["e"][3])
            push!(ids, r["id"])
        end
    end
    if isempty(ids)
        return 0
    end
    w = zeros(length(ids))
    for i = eachindex(w)
        x1 = s["receiver"][ids[i]]["n"][1]
        x2 = s["receiver"][ids[i]]["n"][3]
        y1 = s["receiver"][ids[i]]["e"][1]
        y2 = s["receiver"][ids[i]]["e"][3]
        w[i] = min(abs(n-x1), abs(n-x2))^2 + min(abs(e-y1), abs(e-y2))^2
    end
    tid = findmax(w)
    return ids[tid[2]]
end

readencodestation(path::AbstractString) = map(readlines(path)) do x
    tl = split(x)
    tag = String(tl[1])
    lat = parse(Float64, tl[2])
    lon = parse(Float64, tl[3])
    # el = parse(Float64, tl[4])
    sname = String.(tl[4:end])
    (tag, lat, lon, sname)
end

function parseymd(s::String)
    if length(s) >= 4
        year = parse(Int, s[1:4])
    else
        return DateTime(2000)
    end
    if length(s) >= 6
        mth = parse(Int, s[5:6])
    else
        return DateTime(year)
    end
    if length(s) >= 8
        d = parse(Int, s[7:8])
    else
        return DateTime(year, mth)
    end
    if length(s) >= 10
        h = parse(Int, s[9:10])
    else
        return DateTime(year, mth, d)
    end
    if length(s) >= 12
        mi = parse(Int, s[11:12])
    else
        return DateTime(year, mth, d, h)
    end
    if length(s) >= 14
        sec = parse(Int, s[13:14])
    else
        return DateTime(year, mth, d, h, mi)
    end
    if length(s) >= 15
        msec = parse(Int, s[15:15])
    else
        return DateTime(year, mth, d, h, mi, sec)
    end
    return DateTime(year, mth, d, h, mi, sec, msec*100)
end