using ArgumentProcessor
addopt!("from", abbr="S", fmt=" %s", required=true)
addopt!("to", abbr="T", fmt=" %s", required=true)
addopt!("pattern", abbr="P", fmt=" %s", required=true)
addopt!("run", abbr="R", fmt=" %s", default=" NONE")
addopt!("nproc", abbr="N", fmt=" %d", required=true)
addflag!("resolution"; abbr="H", help="true is high resolution")
IN = ArgumentProcessor.parse(ARGS)

srcdir = abspath(IN.from)
tgtdir = abspath(IN.to)
if !isdir(tgtdir)
    mkpath(tgtdir)
end
midfix = IN.pattern
np = IN.nproc
if IN.run == "NONE"
    t = splitpath(srcdir)
    idx = findfirst(==("OUTPUT_FILES"), t)
    rundir = joinpath(t[1:idx-1])
else
    rundir = abspath(IN.run)
end

@info "$midfix in $srcdir -> $tgtdir\nrunning in $rundir"

files = filter(x -> startswith(x, "proc") &&
                    endswith(x, ".bin") &&
                    contains(x, midfix),
               readdir(srcdir))
if length(files) < 1
    @error "No file exist"
    exit(0)
end
its = (files .|> x -> x[end-9:end-4]) |> sort |> unique
Threads.@threads for it in its
    cmd = Cmd(Cmd(["./bin/xcombine_vol_data_vtk", "0", string(np - 1),
                    midfix * "_it" * it, srcdir, tgtdir, IN.resolution ? "1" : "0"]); dir = rundir)
    run(pipeline(cmd, devnull))
end
