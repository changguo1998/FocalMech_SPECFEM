include("env.jl")
tgt = abspath(PTH["batch2"])
frm = abspath(PTH["raw2"])

evts = readdir(frm)

for e = evts
    epath = joinpath(tgt, String(e[1:end-4]))
    if !isdir(epath)
        mkpath(epath)
    end
    cmd = Cmd(`Event2SAC.sh -b $(joinpath(frm, e)) $epath`)
    try
        run(cmd)
    catch
        println(' ')
    end
end