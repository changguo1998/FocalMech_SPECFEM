using Distributed, ArgumentProcessor, TOML, Mmap, Printf

const INSLURM = "SLURM_JOBID" in keys(ENV)

if INSLURM
    using ClusterManagers
end

addopt!("settingfile"; abbr="S", default=" " * joinpath(pwd(), "setting.toml"),
        fmt=" %s", help="path to setting file")
addopt!("nproc"; abbr="P", default=" 1", fmt=" %d", help="number of processor")
addopt!("nthread"; abbr="T", default=" 1", fmt=" %d", help="number of thread")

input = ArgumentProcessor.parse(ARGS)

if input.nthread > 1
    @info "add $(input.nproc) processors with threading $(input.nthread)"
    addprocs(input.nproc; exeflags="-t " * string(input.nthread))
else
    addprocs(input.nproc)
end

@everywhere begin
    include("../include/io.jl")
    include("../include/SEMshell.jl")
    function calproc(pars::Tuple{String,String,String}, rchannel::RemoteChannel)
        gn = pars[1]
        lhead = pars[2]
        lfield = pars[3]
        (_, _, LREC, _, lglob2glob, lrec2rec, ξr, ηr, γr) = open(SEMshell.readlocalhead, lhead)
        if isempty(lrec2rec)
            sleep(round(rand(), digits=2))
            put!(rchannel, myid())
            return (Int[], zeros(Float32, 1, 6, 0))
        end
        (_, _, NSPEC, _, NGLL, NREC, spec2glob, rec2spec, l, network, station, hp, νr) = open(SEMshell.readglob, gn)
        (NDIM, NT, LGLOB, lveloc) = open(SEMshell.readlocalfield, lfield)
        loop_on_rec = NSPEC > NREC
        g = zeros(Float32, NT, 6, LREC)
        dvidxj = zeros(Float32, NDIM, NDIM, LREC)
        id = zeros(Int, LREC)
        Threads.@threads for irl = 1:LREC
            ir = lrec2rec[irl]
            id[irl] = SEMshell.decodestation(network[ir] * station[ir])
        end
        if loop_on_rec
            (tmp, tmp1, tmp2, tmp3, glob2lglob) = SEMshell.prepare_∇rec(NDIM, NGLL, LREC, LGLOB, lglob2glob)
        else
            (tmp, tmp1, tmp2, tmp3, LSPEC, glob2lglob, lspec2spec, spec2lspec) = SEMshell.prepare_∇spec(NDIM, NSPEC,
                                                                                                        NGLL, LREC,
                                                                                                        lrec2rec,
                                                                                                        LGLOB,
                                                                                                        lglob2glob,
                                                                                                        NREC,
                                                                                                        rec2spec)
        end
        println(lfield, " begin")
        for it = 1:NT
            for idim = 1:NDIM
                if loop_on_rec
                    SEMshell.∇andinterp_looprec!(@view(dvidxj[idim, :, :]), @view(lveloc[idim, :, it]), NDIM,
                                                 NSPEC, NGLL, NREC, rec2spec, spec2glob, νr, hp, l, LREC, LGLOB,
                                                 lrec2rec, ξr, ηr, γr, glob2lglob, tmp, tmp1, tmp2, tmp3)
                else
                    SEMshell.∇andinterp_loopspec!(@view(dvidxj[idim, :, :]), @view(lveloc[idim, :, it]), NDIM,
                                                  NSPEC, NGLL, NREC, rec2spec, spec2glob, νr, hp, l, LREC, LGLOB,
                                                  lrec2rec, ξr, ηr, γr, glob2lglob, spec2lspec, LSPEC, lspec2spec,
                                                  tmp, tmp1, tmp2, tmp3)
                end
            end
            Threads.@threads for irl = 1:LREC
                g[it, 1, irl] = dvidxj[2, 2, irl] # 11
                g[it, 2, irl] = dvidxj[1, 1, irl] # 22
                g[it, 3, irl] = dvidxj[3, 3, irl] # 33
                g[it, 4, irl] = (dvidxj[1, 2, irl] + dvidxj[2, 1, irl]) # 12
                g[it, 5, irl] = -(dvidxj[3, 2, irl] + dvidxj[2, 3, irl]) # 13
                g[it, 6, irl] = -(dvidxj[3, 1, irl] + dvidxj[1, 3, irl]) # 23
            end
        end
        println(lfield, " end")
        put!(rchannel, myid())
        return (id, g)
    end
end

function main()
    println("Load setting")
    shotsetting = SEMshell.loadshotsetting(abspath(input.settingfile))
    stationInfoPath = abspath(shotsetting[4].wkdir, "OUTPUT_FILES")
    Nsemproc = shotsetting[4].nproc_ξ * shotsetting[4].nproc_η
    gfile = normpath(stationInfoPath, "globinterp.bin")

    joblist = Vector{Tuple{String,String,String}}(undef, 3 * Nsemproc)
    jobcmp = zeros(Int, 3 * Nsemproc)
    signchannel = RemoteChannel(() -> Channel{Int}(input.nproc))
    for isetting = 1:3
        if shotsetting[isetting].forceDirc == 'N'
            icmp = 1
        elseif shotsetting[isetting].forceDirc == 'E'
            icmp = 2
        elseif shotsetting[isetting].forceDirc == 'D'
            icmp = 3
        end
        for i = 1:Nsemproc
            joblist[(isetting-1)*Nsemproc+i] = (gfile,
                                                abspath(stationInfoPath, @sprintf("proc%06dlocalinterp.bin", i - 1)),
                                                abspath(shotsetting[isetting].wkdir, "OUTPUT_FILES",
                                                        @sprintf("proc%06ddispl.bin", i - 1)))
            jobcmp[(isetting-1)*Nsemproc+i] = icmp
        end
    end

    println("Submit job")
    jobid = zeros(Int, input.nproc)
    jobbuffer = Vector{Future}(undef, length(joblist))
    jobfinish = falses(length(joblist))
    proc_avail = Channel{Int}(input.nproc)
    lk = ReentrantLock()

    submitjob = @task begin
        currentjob = 1
        while currentjob <= length(joblist)
            iproc = take!(proc_avail)
            lock(lk)
            try
                jobbuffer[iproc-1] = remotecall(calproc, iproc, joblist[currentjob], signchannel)
                jobid[iproc-1] = currentjob
                println("  submit job ", currentjob, " to worker ", iproc)
                currentjob += 1
            finally
                unlock(lk)
            end
        end
        println("All jobs are submitted")
    end
    println("Check worker response")
    for i = 1:input.nproc
        println("  worker ", i + 1, " resp: ", remotecall_fetch(myid, i + 1))
    end
    schedule(submitjob)
    for i = 1:input.nproc
        put!(proc_avail, i + 1)
    end

    println("Calculate coor")
    # coorN = (shotsetting[4].receiver_n[1]:shotsetting[4].receiver_n[2]:shotsetting[4].receiver_n[3]) .-
    #         shotsetting[4].modelrange[1] / 2.0
    # coorE = (shotsetting[4].receiver_e[1]:shotsetting[4].receiver_e[2]:shotsetting[4].receiver_e[3]) .-
    #         shotsetting[4].modelrange[2] / 2.0
    # coorD = (shotsetting[4].receiver_d[1]:shotsetting[4].receiver_d[2]:shotsetting[4].receiver_d[3])
    # Ln = length(coorN)
    # Le = length(coorE)
    # Ld = length(coorD)
    # Lt = shotsetting[1].outnpts
    # idtable = zeros(Int, 3, Ld * Le * Ln)
    # for ix = 1:Le, iy = 1:Ln, iz = 1:Ld
    #     ir = SEMshell.stationnumber(ix, iy, iz, Le, Ln, Ld)
    #     idtable[1, ir] = iz
    #     idtable[2, ir] = ix
    #     idtable[3, ir] = iy
    # end
    rg = shotsetting[4].receiver
    stable = SEMshell.stationtable(rg)

    println("Collect data")
    # io = open(abspath(shotsetting[1].wkdir, "..", "glib.bin"), "w+")
    # write(io, Float32(shotsetting[1].risetime))
    # write(io, Int32(Ln))
    # write(io, Int32(Le))
    # write(io, Int32(Ld))
    # write(io, Int32(Lt))
    # write(io, Float32.(coorN))
    # write(io, Float32.(coorE))
    # write(io, Float32.(coorD))
    # write(io, Float32.(range(0.0; step=shotsetting[1].outdt, length=Lt)))
    iovec = Vector{IO}(undef, length(rg))
    Gvec = Vector{Array{Float32,6}}(undef, length(rg))
    igtable = Dict{Int,Int}()
    Lt = shotsetting[1].outnpts
    for i in eachindex(rg)
        igtable[rg[i].id] = i
        iovec[i] = open(abspath(shotsetting[1].wkdir, "..", @sprintf("glib_tmp1_%d.bin", rg[i].id)), "w+")
        cN = (rg[i].n[1]:rg[i].n[2]:rg[i].n[3]) .- shotsetting[1].srcloc[1]
        cE = (rg[i].e[1]:rg[i].e[2]:rg[i].e[3]) .- shotsetting[1].srcloc[2]
        cD = rg[i].d[1]:rg[i].d[2]:rg[i].d[3]
        write(iovec[i], Float32(shotsetting[1].risetime))
        write(iovec[i], Int32(length(cN)))
        write(iovec[i], Int32(length(cE)))
        write(iovec[i], Int32(length(cD)))
        write(iovec[i], Int32(Lt))
        write(iovec[i], Float32.(cN))
        write(iovec[i], Float32.(cE))
        write(iovec[i], Float32.(cD))
        write(iovec[i], Float32.(range(0.0; step=shotsetting[1].outdt, length=Lt)))
        Gvec[i] = Mmap.mmap(iovec[i], Array{Float32,6}, (Lt, 6, 3, length(cD), length(cE), length(cN)))
    end

    # G = Mmap.mmap(io, Array{Float32,6}, (Lt, 6, 3, Ld, Le, Ln))
    while !all(jobfinish)
        println("  wait for data")
        iproc = take!(signchannel)
        jid = 0
        ids = Int[]
        g = zeros(Float32,0,0,0)
        lock(lk)
        try
            jid = jobid[iproc-1]
            (ids, g) = fetch(jobbuffer[iproc-1])
            println("Get data for job ", jid, " from worker ", iproc)
            put!(proc_avail, iproc)
        finally
            unlock(lk)
        end
        if isempty(ids)
            jobfinish[jid] = true
            continue
        end
        icmp = jobcmp[jid]
        # Threads.@threads for irl in eachindex(ids)
        #     id = ids[irl]
        #     iz = idtable[1, id]
        #     iy = idtable[2, id]
        #     ix = idtable[3, id]
        #     G[:, :, icmp, iz, iy, ix] .= g[:, :, irl]
        # end
        Threads.@threads for irl in eachindex(ids)
            id = ids[irl]
            ig = igtable[stable[id, 2]]
            iz = stable[id, 3]
            iy = stable[id, 4]
            ix = stable[id, 5]
            Gvec[ig][:, :, icmp, iz, iy, ix] .= g[:, :, irl]
        end
        jobfinish[jid] = true
        println("job ", jid, " finish")
    end
    # Mmap.sync!(G)
    # close(io)
    for i in eachindex(iovec)
        Mmap.sync!(Gvec[i])
        close(iovec[i])
    end
    return nothing
end

main()