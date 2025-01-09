#!/usr/bin/env julia
using ArgumentProcessor, Dates

function getetime(path::AbstractString)
    if isfile(path)
        lines = readlines(path)
        idx1 = findlast(contains('%'), lines)
        t = replace(replace(lines[idx1], "We have done" => ""), "% of that" => "")
        s = parse(Float64, strip(t))
        if s < 100.0
            idx2 = findlast(contains("The run will finish approximately"), lines)
            f = replace(lines[idx2], "The run will finish approximately on (in local time):" => "")
        else
            f = "--:--"
        end
        return (s, strip(f))
    else
        return (0.0, "err")
    end
end

function getlaststamp(path::AbstractString)
    if isdir(path)
        fs = filter(contains("timestamp"), readdir(path))
        it = map(fs) do f
            parse(Int, strip(replace(f, "timestamp" => "")))
        end
        if isempty(it)
            return ""
        end
        (_, id) = findmax(it)
        return joinpath(path, fs[id])
    else
        return ""
    end
end

addopt!("wkdir"; abbr="D", fmt=" %s", required=true)
addopt!("waittime"; abbr="N", fmt=" %d", required=true)
input = ArgumentProcessor.parse(ARGS)
dirs = ("N", "E", "D")

while true
    stat = map(dirs) do d
        fp = getlaststamp(joinpath(input.wkdir, "shot"*d, "OUTPUT_FILES"))
        t = getetime(fp)
        (d, t[1], t[2])
    end
    println("\n", now())
    for (dirc, pernt, et) in stat
        println(dirc, ": ", pernt, "%, finish at ", et)
    end
    if all(map(v -> v[2] == 100.0, stat))
        break
    end
    sleep(input.waittime)
end
