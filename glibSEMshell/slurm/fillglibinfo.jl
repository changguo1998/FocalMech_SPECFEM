using ArgumentProcessor, TOML, Printf

NLLOCrootpath = abspath(@__DIR__, "../NLLOCshell/")

include(joinpath(NLLOCrootpath, "NLLOCshell.jl"))


addopt!("settingfile"; abbr="S", fmt=" %s", default=" setting.toml")

input = ArgumentProcessor.parse(ARGS)

S = TOML.parsefile(input.settingfile)

# source time function

stf_shot = let
    t = parse.(Float64, readlines(joinpath(S["wkdir"], "shotD/DATA/EXTERNAL_SOURCE_TIME_FUNCTION")))
    dt = t[1]
    Float32.(t[2:end])
end
Lstfshot = length(stf_shot)
stf = zeros(Float32, S["npts"])
for it = eachindex(stf)
    stf[it] = stf_shot[round(Int, (it-1)/(S["npts"]-1)*(Lstfshot-1))+1]
end

cmd = Cmd(Cmd(["julia",
    joinpath(NLLOCrootpath, "calculateTravelTime.jl"),
    "-M", S["modelfile"],
    "-C", joinpath(S["wkdir"], "tlibvar"),
    "-T", joinpath(S["wkdir"], "tlib.bin"),
    "-S", input.settingfile,
    "-H", string(S["tlib_h"])
]); dir=S["wkdir"])

if !isfile(joinpath(S["wkdir"], "tlib.bin"))
    run(cmd)
end

(dxt, dyt, dzt, oxt, oyt, ozt, tt) = open(NLLOCshell.ttlib_readall, joinpath(S["wkdir"], "tlib.bin"))
nxt = size(tt, 4)
nyt = size(tt, 3)
nzt = size(tt, 2)

for iglib = eachindex(S["receiver"])
    rg = S["receiver"][iglib]
    glibid = rg["id"]
    glibpath = joinpath(S["wkdir"], @sprintf("glib_%d.bin", glibid))
    if !isfile(glibpath)
        continue
    end
    io = open(glibpath, "r")
    flag = read(io, UInt8)
    partstart = zeros(UInt64, 4); read!(io, partstart)
    seek(io, 1+8*4+4*3)
    nt = read(io, Int32)
    dt = read(io, Float32)
    seek(io, partstart[1])
    ns = zeros(Int32, 3); read!(io, ns)
    xs = zeros(Float32, 6); read!(io, xs)
    (nx, ny, nz) = ns
    (x0, y0, z0, dx, dy, dz) = xs
    tracepospos = read(io, UInt64)
    seek(io, tracepospos)
    grid = zeros(UInt64, nz, ny, nx); read!(io, grid);
    close(io)
    io = open(glibpath, "a")
    seek(io, 1+8*4)
    write(io, Float32(S["srclat"]))
    write(io, Float32(S["srclon"]))
    write(io, Float32(S["srcdep"]))
    flush(io)
    seek(io, 1+8*4+4*3+4+4)
    write(io, stf)
    flush(io)

    for ix = 1:nx, iy = 1:ny, iz = 1:nz
        x = x0 + (ix-1)*dx
        y = y0 + (iy-1)*dy
        z = z0 + (iz-1)*dz
        ixt = floor(Int, (x - oxt) / dxt) + 1
        iyt = floor(Int, (y - oyt) / dyt) + 1
        izt = floor(Int, (z - ozt) / dzt) + 1
        if (ixt < 1) || (ixt > nxt)
            error("Out of index range")
        end
        if (iyt < 1) || (iyt > nyt)
            error("Out of index range")
        end
        if (izt < 1) || (izt > nzt)
            error("Out of index range")
        end
        if ixt == nxt
            ixt = nxt - 1
        end
        if iyt == nyt
            iyt = nyt - 1
        end
        if izt == nzt
            izt = nzt - 1
        end
        hxt = (x - oxt) / dxt - ixt + 1;
        hyt = (y - oyt) / dyt - iyt + 1;
        hzt = (z - ozt) / dzt - izt + 1;
        tp = tt[1,izt,   iyt,   ixt]   * (1.0 - hxt) * (1.0 - hyt) * (1.0 - hzt) +
            tt[1, izt+1, iyt,   ixt]   * (1.0 - hxt) * (1.0 - hyt) * hzt         +
            tt[1, izt,   iyt+1, ixt]   * (1.0 - hxt) * hyt         * (1.0 - hzt) +
            tt[1, izt+1, iyt+1, ixt]   * (1.0 - hxt) * hyt         * hzt         +
            tt[1, izt,   iyt,   ixt+1] * hxt         * (1.0 - hyt) * (1.0 - hzt) +
            tt[1, izt+1, iyt,   ixt+1] * hxt         * (1.0 - hyt) * hzt         +
            tt[1, izt,   iyt+1, ixt+1] * hxt         * hyt         * (1.0 - hzt) +
            tt[1, izt+1, iyt+1, ixt+1] * hxt         * hyt         * hzt;
        ts =tt[2,izt,   iyt,   ixt]   * (1.0 - hxt) * (1.0 - hyt) * (1.0 - hzt) +
            tt[2, izt+1, iyt,   ixt]   * (1.0 - hxt) * (1.0 - hyt) * hzt         +
            tt[2, izt,   iyt+1, ixt]   * (1.0 - hxt) * hyt         * (1.0 - hzt) +
            tt[2, izt+1, iyt+1, ixt]   * (1.0 - hxt) * hyt         * hzt         +
            tt[2, izt,   iyt,   ixt+1] * hxt         * (1.0 - hyt) * (1.0 - hzt) +
            tt[2, izt+1, iyt,   ixt+1] * hxt         * (1.0 - hyt) * hzt         +
            tt[2, izt,   iyt+1, ixt+1] * hxt         * hyt         * (1.0 - hzt) +
            tt[2, izt+1, iyt+1, ixt+1] * hxt         * hyt         * hzt;
        seek(io, grid[iz,iy,ix])
        write(io, Float32(tp))
        write(io, Float32(ts))
        flush(io)
    end

    close(io)
end