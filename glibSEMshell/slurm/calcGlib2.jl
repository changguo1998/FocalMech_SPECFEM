using ArgumentProcessor, TOML, Printf

include("../include/io.jl")
include("../include/SEMshell.jl")

addopt!("settingfile"; abbr="S", default=" " * joinpath(pwd(), "setting.toml"),
        fmt=" %s", help="path to setting file")
addopt!("nthread"; abbr="T", default=" 1", fmt=" %d", help="number of thread")

input = ArgumentProcessor.parse(ARGS)
shotsetting = SEMshell.loadshotsetting(abspath(input.settingfile))
rootPath = abspath(shotsetting[4].wkdir, "..")
Nsemproc = shotsetting[4].nproc_ξ * shotsetting[4].nproc_η
rg = shotsetting[4].receiver
stable = SEMshell.stationtable(rg)
open(joinpath(rootPath, "station_id_table.bin"), "w") do io
    write(io, Int32(size(stable, 1)))
    write(io, Int32.(stable));
end    

Lt = shotsetting[1].outnpts
for i in eachindex(rg)
    io = open(abspath(rootPath, @sprintf("glib_tmp2_%d.bin", rg[i].id)), "w")
    cN = (rg[i].n[1]:rg[i].n[2]:rg[i].n[3]) .- shotsetting[1].srcloc[1]
    cE = (rg[i].e[1]:rg[i].e[2]:rg[i].e[3]) .- shotsetting[1].srcloc[2]
    cD = rg[i].d[1]:rg[i].d[2]:rg[i].d[3]
    write(io, Float32(shotsetting[1].risetime))
    write(io, Int32(length(cN)))
    write(io, Int32(length(cE)))
    write(io, Int32(length(cD)))
    write(io, Int32(Lt))
    write(io, Float32.(cN))
    write(io, Float32.(cE))
    write(io, Float32.(cD))
    write(io, Float32.(range(0.0; step=shotsetting[1].outdt, length=Lt)))
    close(io);
end

glib_calc_path = abspath(SEMshell.SEMHOME, "bin/calcglib")
cmd = Cmd([glib_calc_path, rootPath, string(Nsemproc), string(length(rg)), string(input.nthread)]);
println(cmd)
run(Cmd(cmd; dir=pwd()))


for i in eachindex(rg)
    ghdrpath = joinpath(rootPath, @sprintf("glib_tmp2_%d.bin", i));
    gdatpath = joinpath(rootPath, @sprintf("glib_tmp3_%d.bin", i));
    gpath = joinpath(rootPath, @sprintf("glib_tmp1_%d.bin", i));
    run(pipeline(Cmd(Cmd(["cat", ghdrpath, gdatpath]); dir=rootPath), gpath));
    rm(ghdrpath)
    rm(gdatpath)
end
