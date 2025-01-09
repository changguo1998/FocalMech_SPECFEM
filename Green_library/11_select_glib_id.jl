using TOML, SeisTools, Printf

centerlon = 102.25
centerlat = 29.75

glibdir = abspath("../../dat/glib_hpc/")
glibs = filter(s->isfile(joinpath(glibdir, s, "setting.toml")), 
    setdiff(readdir(glibdir), ["generate_index.jl", "index.bak.txt", "index.txt"]))

glibid = map(glibs) do g
    s = TOML.parsefile(joinpath(glibdir, g, "setting.toml"))
    dist = SeisTools.Geodesy.distance(s["srclat"], s["srclon"], centerlat, centerlon) * 0.001
    baz = SeisTools.Geodesy.azimuth(s["srclat"], s["srclon"], centerlat, centerlon)
    n = dist*cosd(baz)
    e = dist*sind(baz)
    x = n + s["srcloc"][1]
    y = e + s["srcloc"][2]
    rid = findfirst(s["receiver"]) do r
        (x > r["n"][1]) &&
        (x < r["n"][3]) &&
        (y > r["e"][1]) &&
        (y < r["e"][3])
    end
    s["receiver"][rid]["id"]
end

for i = eachindex(glibs)
    @printf("%s glib_%d.bin\n", glibs[i], glibid[i])
end

mkpath("temporary_glib")
for i = eachindex(glibs)
    gdir = joinpath("temporary_glib", glibs[i])
    mkpath(gdir)
    gf = @sprintf("glib_%d.bin", glibid[i])
    if !isfile(joinpath(gdir, gf))
        cp(joinpath(glibdir, glibs[i], gf), joinpath(gdir, gf))
    end
    if !isfile(joinpath(gdir, "setting.toml"))
        cp(joinpath(glibdir, glibs[i], "setting.toml"), joinpath(gdir, "setting.toml"))
    end
end
