include("env.jl")

batch1dir = PTH["dat_selectstation"]
es = readdir(batch1dir)
respallfile = joinpath(PTH["root"], "dat", "RESP.all")

Threads.@threads for e in es
    open(joinpath(batch1dir, e, "sac", "rmresp.bash"), "w") do io
        println(io, "sac <<EOF")
        println(io, "r *.SAC")
        println(io, "rmean; rtr; taper;")
        println(io, "trans from evalresp fname $respallfile to none freq 0.01 0.02 30 40")
        println(io, "w over")
        println(io, "q")
        println(io, "EOF")
    end
    run(
        pipeline(
            Cmd(Cmd(["bash", "rmresp.bash"]); dir=joinpath(batch1dir, e, "sac")),
            devnull
            )
    )
    rm(joinpath(batch1dir, e, "sac", "rmresp.bash"))
end