using Printf, Downloads

function srtm_download_block(lat::Real, lon::Real)
    if (lat <= -60) || (lat >= 60)
        error("lat out of bound")
    end
    if (lon <= -180) || (lon >= 180)
        error("lon out of bound")
    end
    Nolat = floor(Int, (60 - lat) / 5.0) + 1
    Nolon = floor(Int, (180 + lon) / 5.0) + 1
    blockfilename = @sprintf("srtm_%02d_%02d.zip", Nolon, Nolat)
    url_s = "https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/ASCII/" * blockfilename
    dat_s = joinpath(@__DIR__, "raw_data", blockfilename)
    if !isfile(dat_s)
        @info "downloading block data $blockfilename"
        mkpath(joinpath(@__DIR__, "raw_data"))
        try
            Downloads.download(url_s, dat_s)
        catch
            @error "Failed to download $blockfilename"
        end
    end
    return nothing
end

if length(ARGS) < 4
    println("Usage:")
    println("    julia download_srtm.jl minlat maxlat minlon maxlon")
    exit(0)
end

minlat = parse(Float64, ARGS[1])
maxlat = parse(Float64, ARGS[2])
minlon = parse(Float64, ARGS[3])
maxlon = parse(Float64, ARGS[4])

lats = range(minlat, stop = maxlat, step = 0.1)
lons = range(minlon, stop = maxlon, step = 0.1)

for lat = lats, lon = lons
    srtm_download_block(lat, lon)
end
