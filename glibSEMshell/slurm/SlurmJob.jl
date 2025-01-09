module SlurmJob

struct ResourceSetting
    jobname::String
    ntask::Int
    nnode::Int
    queue::String
    ntaskpernode::Int
    cpupertask::Int
    time::String
    output::String
    error::String
    mempernode::String
    mempercpu::String
end

# function ResourceSetting(jobname::AbstractString="defaultjob", ntask::Integer=0, nnode::Integer=0,
#                        queue::AbstractString="hfacnormal01", ntaskpernode::Integer=0, cpupertask::Integer=0,
#                        time::AbstractString="", output::AbstractString="%j.out", error::AbstractString="%j.err",
#                        mempernode::AbstractString="", mempercpu::AbstractString="")
#     return ResourceSetting(String(jobname), Int(ntask), Int(nnode), String(queue), Int(ntaskpernode), Int(cpupertask),
#                          String(time), String(output), String(error), String(mempernode), String(mempercpu))
# end

function ResourceSetting(; jobname::AbstractString="defaultjob", ntask::Integer=0, nnode::Integer=0,
                         queue::AbstractString="hfacnormal01", ntaskpernode::Integer=0, cpupertask::Integer=0,
                         time::AbstractString="", output::AbstractString="%j.out", error::AbstractString="%j.err",
                         mempernode::AbstractString="", mempercpu::AbstractString="")
    # if iszero(ntaskpernode)
    #     ntaskpernode = ceil(Int, ntask/nnode)
    # end
    return ResourceSetting(String(jobname), Int(ntask), Int(nnode), String(queue), Int(ntaskpernode), Int(cpupertask),
                           String(time), String(output), String(error), String(mempernode), String(mempercpu))
end

function writescript(path::AbstractString, setting::ResourceSetting, cmd::Vector{<:AbstractString})
    # if setting.ntask != (setting.ntaskpernode*setting.nnode)
    #     error("ntask not equal to nnode*ntaskpernode")
    # end
    open(path, "w") do io
        println(io, "#!/bin/bash")
        println(io, "#SBATCH --job-name=", setting.jobname)
        if !iszero(setting.ntask)
            println(io, "#SBATCH --ntasks=", setting.ntask)
        end
        if !iszero(setting.nnode)
            println(io, "#SBATCH --nodes=", setting.nnode)
        end
        println(io, "#SBATCH --partition=", setting.queue)
        if !iszero(setting.ntaskpernode)
            println(io, "#SBATCH --ntasks-per-node=", setting.ntaskpernode)
        end
        if !iszero(setting.cpupertask)
            println(io, "#SBATCH --cpus-per-task=", setting.cpupertask)
        end
        if !isempty(setting.time)
            println(io, "#SBATCH --time=", setting.time)
        end
        if !isempty(setting.output)
            println(io, "#SBATCH --output=", setting.output)
        end
        if !isempty(setting.error)
            println(io, "#SBATCH --error=", setting.error)
        end
        if !isempty(setting.mempernode)
            println(io, "#SBATCH --mem=", setting.mempernode)
        end
        if !isempty(setting.mempercpu)
            println(io, "#SBATCH --mem-per-cpu=", setting.mempercpu)
        end
        for l in cmd
            println(io, l)
        end
    end
    return nothing
end

function writescript(path::AbstractString, setting::ResourceSetting, cmd::AbstractString)
    writescript(path, setting, [cmd])
end

end
