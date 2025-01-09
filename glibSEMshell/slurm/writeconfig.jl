#!/usr/bin/env julia
using Dates, Printf, TOML, ArgumentProcessor

include(abspath(@__DIR__, "../include/SEMshell.jl"))
include(abspath(@__DIR__, "SlurmJob.jl"))

function main()
    addopt!("settingfile"; abbr="S", default=" " * joinpath(pwd(), "setting.toml"),
            fmt=" %s", help="path to setting file")
    addopt!("nodes"; abbr="N", default=" 1", fmt=" %d", help="Number of nodes")
    addopt!("label"; abbr="L", default=" default", fmt=" %s")
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
        SEMshell.writestationposition(shotsetting[i]); GC.gc()
        SEMshell.writesourcetimefunction(shotsetting[i])
        SEMshell.writemainpar(shotsetting[i]); GC.gc()
        SEMshell.writeforcesolution(shotsetting[i])
    end

    for i = [2, 3]
        cp(shotsetting[1].wkdir, shotsetting[i].wkdir)
    end

    for s in shotsetting
        nproc = s.nproc_ξ * s.nproc_η
        SEMshell.writeforcesolution(s) # ! cannot be removed
        hpcsetting = SlurmJob.ResourceSetting(;
                                              jobname="shot_" * input.label * "_" * s.forceDirc,
                                              ntask=nproc,
                                              nnode=input.nodes)
        SlurmJob.writescript(joinpath(s.wkdir, "submitjob.slurm"), hpcsetting,
                             ["mpirun -np " * string(nproc) * " ./bin/xmeshfem3D && \\",
                              "mpirun -np " * string(nproc) * " ./bin/xgenerate_databases && \\",
                              "mpirun -np " * string(nproc) * " ./bin/xspecfem3D && \\",
                              "touch ../shot"*s.forceDirc*".flag"])
    end

    @info "Done at $(now())"
    return nothing
end

main()
