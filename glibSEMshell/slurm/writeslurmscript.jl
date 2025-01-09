#!/usr/bin/env julia
using ArgumentProcessor, Printf

include("SlurmJob.jl")

grp1 = Group("config";
             opts=Option[Option("configthreads"; abbr="Tc", fmt=" %d", default=" 32",
                                help="Threads when writing configs"),
                         Option("forwardnodes"; abbr="Ns", fmt=" %d", default=" 1",
                                help="Nodes when running sem")])
grp2 = Group("glib";
             opts=Option[Option("glibnodes"; abbr="Ng", fmt=" %d", default=" 1",
                                help="Nodes when calculate green lib"),
                        #  Option("glibprocs"; abbr="Pg", fmt=" %d", default=" 4",
                        #         help="Processors when calculate green lib"),
                         Option("glibthreads"; abbr="Tg", fmt=" %d", default=" 16",
                                help="Threads when calculate green lib"),
                         Option("compressthreads"; abbr="Cg", fmt=" %d", default=" 16",
                                help="Threads when compress green lib"),
                         Option("compressExpbit"; abbr="Cb", fmt=" %d", default=" 3",
                                help="expbit"),
                         Option("compressSigbit"; abbr="Cs", fmt=" %d", default=" 7",
                                help="sigbit")])
grp3 = Group("glob";
             opts=Option[Option("settingfile"; abbr="S", fmt=" %s", default=" " * abspath(pwd(), "setting.toml")),
             Option("label"; abbr="L", fmt=" %s", default=" default")])

input = ArgumentProcessor.parse(ARGS, [grp1, grp2, grp3])

@info @sprintf("Configuration file: %s\nWhile writting config, using threads: %d\nRunning SEM on %d node(s)
While calculating glib, using %d thread(s) on %d node(s)\n",
               input.glob.settingfile, input.config.configthreads, input.config.forwardnodes, input.glib.glibthreads, input.glib.glibnodes)

writerpath = abspath(@__DIR__, "writeconfig.jl")
cglibpath = abspath(@__DIR__, "calcGlib2.jl")
compressglibpath = abspath(@__DIR__, "compress.jl")
fillinfoglibpath = abspath(@__DIR__, "fillglibinfo.jl")

SlurmJob.writescript(joinpath(pwd(), "writeconfig.slurm"),
                     SlurmJob.ResourceSetting(; jobname="writeconfig_"*input.glob.label,
                                              nnode=1,
                                              ntask=1,
                                              cpupertask=input.config.configthreads),
                     [join(["julia", "-t", string(input.config.configthreads), writerpath,
                           "-S", abspath(input.glob.settingfile),
                           "-N", string(input.config.forwardnodes),
                           "-L", input.glob.label], ' ') * " && \\",
                    "touch writeconfig.flag"])

SlurmJob.writescript(joinpath(pwd(), "calcglib.slurm"),
                     SlurmJob.ResourceSetting(; jobname="glib_"*input.glob.label,
                                              nnode=input.glib.glibnodes,
                                              ntask=1,
                                              cpupertask=input.glib.glibthreads),
                     [join(["julia", cglibpath,
                           "-S", abspath(input.glob.settingfile),
                           "-T", string(input.glib.glibthreads)], ' ') * " && \\",
                    "touch calcglib.flag"])

SlurmJob.writescript(joinpath(pwd(), "compressglib.slurm"),
    SlurmJob.ResourceSetting(; jobname="cmrs_"*input.glob.label,
        nnode=input.glib.glibnodes,
        ntask=1,
        cpupertask=input.glib.compressthreads),
    [join(["julia", compressglibpath,
        "-S", input.glob.settingfile,
        "-T", string(input.glib.compressthreads),
        "-Be", input.glib.compressExpbit,
        "-Bs", input.glib.compressSigbit], ' ') * " && \\",
    "touch compress.flag"]
)

SlurmJob.writescript(joinpath(pwd(), "fillinfo.slurm"),
    SlurmJob.ResourceSetting(; jobname="fill_info_"*input.glob.label,
        nnode=1,
        ntask=1,
        cpupertask=4),
    [join(["julia", fillinfoglibpath,
        "-S", input.glob.settingfile], ' ') * " && \\",
    "touch fillinfo.flag"]
)

open(joinpath(pwd(), "copydata.sh"), "w") do io
    for c in ['D', 'E', 'N']
        println(io, "cp shotT/OUTPUT_FILES/*specinfo.bin shot", c, "/OUTPUT_FILES/")
    end
end
