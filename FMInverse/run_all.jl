include("env.jl")
batchdir = PTH["dat_selectstation"]

# scriptdir = joinpath(@__DIR__, "inverse_script")

lk = ReentrantLock()
es = readdir(batchdir)
# es = ["201408031630130"]
eflag = falses(length(es))
tasks = Vector{Task}(undef, Threads.nthreads())

for i = 1:Threads.nthreads()
    tasks[i] = @task begin
        while true
            global lk, es, eflag
            lock(lk)
            j = 0
            try
                j = findfirst(!, eflag)
                if !isnothing(j)
                    eflag[j] = true
                end
            finally
                unlock(lk)
            end
            if isnothing(j)
                break
            end
            e = es[j]
            # println(e)
            try
                # run(Cmd(`bash run_event.sh`; dir=srcdir); wait=true)
                if isfile(joinpath(PTH["dat_selectstation"], e, "phase.txt"))
                    run(Cmd(Cmd(["julia", "preprocess.jl", e]); dir=PTH["invscript"]); wait=true)
                    run(Cmd(Cmd(["julia", "-t", "16", "inverse.jl", e]); dir=PTH["invscript"]); wait=true)
                    # run(Cmd(Cmd(["julia", "plot.jl", e]); dir=PTH["invscript"]); wait=true)
                end
            catch x
                lock(lk) do
                    open("batchall_inverse_log.txt", "a+") do io
                        println(io, e)
                    end
                end
            end
            lock(lk) do
                println(e, " finish")
            end
        end
        # println("thread ", Threads.threadid(), " task ", i, " exit")
    end
end

Threads.@threads for t in tasks
    schedule(t)
end

while true
    f = nothing
    lock(lk)
    try
        f = findfirst(t->(!istaskdone(t)) && (!istaskfailed(t)), tasks)
    finally
        unlock(lk)
    end
    if isnothing(f)
        break
    else
        wait(tasks[f])
    end
end
