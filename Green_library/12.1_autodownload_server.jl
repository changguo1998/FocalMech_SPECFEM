function rm_temp_file(gdir::String)
    tmpfiles = String[]
    allfiles = readdir(gdir)
    append!(tmpfiles, filter(startswith("glib_tmp1_"), allfiles))
    append!(tmpfiles, filter(endswith(".slurm"), allfiles))
    append!(tmpfiles, filter(endswith(".out"), allfiles))
    append!(tmpfiles, filter(endswith(".err"), allfiles))
    append!(tmpfiles, filter(endswith(".submited"), allfiles))
    for _d = ["shotD", "shotE", "shotN", "shotT", "tlibvar", "copydata.sh",
        "station_id_table.bin", "tlib.bin"]
        if _d in allfiles
            push!(tmpfiles, _d)
        end
    end

    for f in tmpfiles
        println(joinpath(gdir,f))
        rm(joinpath(gdir,f), recursive=true)
    end
end

glibdir = abspath("../../dat/glib/")
station_can_be_download = String[]
for g in readdir(glibdir)
    if isfile(joinpath(glibdir, g, "fillinfo.flag")) &&
        (!isfile(joinpath(glibdir, g, "downloading.flag"))) &&
        (!isfile(joinpath(glibdir, g, "download.flag")))
        push!(station_can_be_download, g)
    end
    if isfile(joinpath(glibdir, g, "download.flag"))
        for gbinfile in filter(endswith(".bin"), readdir(joinpath(glibdir, g)))
            rm(joinpath(glibdir, g, gbinfile))
        end
    end
end

if !isempty(station_can_be_download)
    println("can be downloaded:")
    foreach(println, station_can_be_download)
end

for g in station_can_be_download
    rm_temp_file(joinpath(glibdir, g))
end

open("station_wait_to_download.txt", "w") do io
    for g in station_can_be_download
        println(io, g)
    end
end
