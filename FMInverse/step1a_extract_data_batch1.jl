include("env.jl")

target = PTH["batch1"]
seedfiles = String[]

for y = readdir(PTH["raw1"])
    if isdir(abspath(PTH["raw1"], y))
        for m = readdir(abspath(PTH["raw1"], y))
            if isdir(abspath(PTH["raw1"], y, m))
                for e = readdir(abspath(PTH["raw1"], y, m))
                    if endswith(e, "SEED")
                        tgtdir = abspath(target, String(e[1:15]))
                        mkpath(tgtdir)
                        cmd = Cmd(["rdseed.linux_64", "-f", abspath(PTH["raw1"], y, m, e), "-d", "-o", "1", "-q", tgtdir])
                        println(cmd)
                        run(cmd)
                        cmd = Cmd(["rdseed.linux_64", "-f", abspath(PTH["raw1"], y, m, e), "-d", "-R", "-q", tgtdir])
                        println(cmd)
                        run(cmd)
                    end
                end
            end
        end
    end
end

