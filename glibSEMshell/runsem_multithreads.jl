#!/usr/bin/env julia
using Dates, Printf, TOML, ArgumentProcessor

include("include/SEMshell.jl")

function main()
    addopt!("settingfile"; abbr="S", default=" " * joinpath(pwd(), "setting.toml"),
            fmt=" %s", help="path to setting file")
    input = ArgumentProcessor.parse(ARGS)
    shotsetting = SEMshell.loadshotsetting(abspath(input.settingfile))

    SEMshell.cleandir(shotsetting[1].wkdir)
    for d in ("DATA", "DATA/meshfem3D_files", "bin", "OUTPUT_FILES", "OUTPUT_FILES/DATABASES_MPI")
        mkpath(joinpath(shotsetting[1].wkdir, d))
    end
    SEMshell.writemeshpar(shotsetting[1])
    SEMshell.genSEMdatabase(shotsetting[1])
    cp(shotsetting[1].wkdir, shotsetting[4].wkdir)

    for i in [1, 4]
        SEMshell.writestationposition(shotsetting[i])
        SEMshell.writesourcetimefunction(shotsetting[i])
        SEMshell.writemainpar(shotsetting[i])
        SEMshell.writeforcesolution(shotsetting[i])
        SEMshell.run_xmeshfem3D(shotsetting[i])
    end

    SEMshell.run_xgenerate_databases(shotsetting[4])
    SEMshell.run_xspecfem3D(shotsetting[4])

    for i = [2, 3]
        cp(shotsetting[1].wkdir, shotsetting[i].wkdir)
    end
    globinfos = filter(readdir(joinpath(shotsetting[4].wkdir, "OUTPUT_FILES"), join=true)) do v
        (_, t) = splitdir(v)
        startswith(t, "proc") && endswith(t, "specinfo.bin")
    end
    for i = 1:3, f in globinfos
        (_, fn) = splitdir(f)
        cp(f, joinpath(shotsetting[i].wkdir, "OUTPUT_FILES", fn))
    end

    Threads.@threads for s in shotsetting[1:3]
        SEMshell.writeforcesolution(s) # ! cannot be removed
        SEMshell.run_xgenerate_databases(s)
        SEMshell.run_xspecfem3D(s)
    end

    @info "Construct Green lib"
    SEMshell.constructgreenlib(shotsetting, abspath(shotsetting[1].wkdir, ".."))
    @info "Done at $(now())"
    return nothing
end

main()
