using Distributed, ArgumentProcessor, TOML

addopt!("datadir"; abbr="C", fmt=" %s", required=true)
addopt!("settingfile"; abbr="S", default=" " * joinpath(pwd(), "setting.toml"),
        fmt=" %s", help="path to setting file")
addopt!("outputdir"; abbr="O", fmt=" %s", required=true, help="Output directory")
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
    using Printf
    # include("include/io.jl")
    # include("include/SEMshell.jl")
    include("/home/guochang/Projects/glibSEMshell/include/io.jl")
    include("/home/guochang/Projects/glibSEMshell/include/SEMshell.jl")
    function calproc(pars::Tuple{String,String,String}, rchannel::RemoteChannel)
        gn = pars[1]
        lhead = pars[2]
        lfield = pars[3]
        (_, _, NSPEC, _, NGLL, NREC, spec2glob, rec2spec, l, network, station, hp, νr) = open(SEMshell.readglob, gn)
        # (_, _, LREC, LGLOB, lglob2glob, lrec2rec, ξr, ηr, γr, lveloc) = open(SEMshell.readlocal, fn)
        (_, _, LREC, _, lglob2glob, lrec2rec, ξr, ηr, γr) = open(SEMshell.readlocalhead, lhead)
        (NDIM, NT, LGLOB, lveloc) = open(SEMshell.readlocalfield, lfield)
        loop_on_rec = NSPEC > NREC
        dvidxj = zeros(Float32, NDIM, NDIM, LREC, NT)
        veloc = zeros(Float32, NDIM, LREC, NT)
        g = zeros(Float32, 6, LREC, NT)
        id = zeros(Int, LREC)
        Threads.@threads for irl = 1:LREC
            ir = lrec2rec[irl]
            id[irl] = SEMshell.decodestation(network[ir] * station[ir])
        end
        (tmpv, glob2lglobv) = SEMshell.prepare_interpfield(NDIM, LREC, LGLOB, lglob2glob)
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
        @info "$(lhead) start"
        for it = 1:NT
            SEMshell.interpolatefield!(@view(veloc[:, :, it]), @view(lveloc[:, :, it]), NDIM, NSPEC, NGLL, NREC,
                                        rec2spec, spec2glob, νr, LREC, LGLOB, lrec2rec, ξr, ηr, γr, glob2lglobv,
                                        tmpv)
            for idim = 1:NDIM
                if loop_on_rec
                    SEMshell.∇andinterp_looprec!(@view(dvidxj[idim, :, :, it]), @view(lveloc[idim, :, it]), NDIM,
                                                    NSPEC, NGLL, NREC, rec2spec, spec2glob, νr, hp, l, LREC, LGLOB,
                                                    lrec2rec, ξr, ηr, γr,
                                                    glob2lglob, tmp, tmp1, tmp2, tmp3)
                else
                    SEMshell.∇andinterp_loopspec!(@view(dvidxj[idim, :, :, it]), @view(lveloc[idim, :, it]), NDIM,
                                                    NSPEC, NGLL, NREC, rec2spec, spec2glob, νr, hp, l, LREC, LGLOB,
                                                    lrec2rec, ξr, ηr, γr,
                                                    glob2lglob, spec2lspec, LSPEC, lspec2spec, tmp, tmp1, tmp2, tmp3)
                end
            end
            Threads.@threads for irl = 1:LREC
                g[1, irl, it] = dvidxj[2, 2, irl, it] # 11
                g[2, irl, it] = dvidxj[1, 1, irl, it] # 22
                g[3, irl, it] = dvidxj[3, 3, irl, it] # 33
                g[4, irl, it] = (dvidxj[1, 2, irl, it] + dvidxj[2, 1, irl, it]) # 12
                g[5, irl, it] = -(dvidxj[3, 2, irl, it] + dvidxj[2, 3, irl, it]) # 13
                g[6, irl, it] = -(dvidxj[3, 1, irl, it] + dvidxj[1, 3, irl, it]) # 23
            end
        end
        @info "$(lhead) end"
        put!(rchannel, myid())
        return (id, veloc, dvidxj, g)
    end
end

function main()
    println("Load setting")
    shotsetting = SEMshell.loadshotsetting(abspath(input.settingfile))
    stationInfoPath = abspath(input.datadir, "../../shotT/OUTPUT_FILES")
    bufFilePath = abspath(input.datadir)
    gfile = normpath(stationInfoPath, "globinterp.bin")
    pfiles = filter(readdir(bufFilePath; join=true)) do v
        (_, t) = splitdir(v)
        startswith(t, "proc") && endswith(t, "displ.bin")
    end
    Nsemproc = length(pfiles)

    joblist = Tuple[]
    signchannel = RemoteChannel(() -> Channel{Int}(input.nproc))
    for i = 1:Nsemproc
        push!(joblist,
              (gfile,
               abspath(stationInfoPath, @sprintf("proc%06dlocalinterp.bin", i - 1)),
               abspath(bufFilePath, @sprintf("proc%06ddispl.bin", i - 1))))
    end
    # for f in pfiles
    #     push!(joblist, (gfile, f))
    # end

    # jobid = eachindex(joblist)
    # jobworker = zeros(Int, length(jobid))
    # for i in jobid
    #     jobworker[i] = mod(i - 1, input.nproc) + 2
    # end
    # resultfuture = Vector{Future}(undef, length(joblist))
    # for i in eachindex(joblist)
    #     resultfuture[i] = remotecall(calproc, mod(i - 1, input.nproc) + 2, joblist[i][1], joblist[i][2])
    # end
    # for i = 1:input.nproc
    #     resultfuture[i] = remotecall(calproc, jobworker[i], joblist[i])
    # end
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
    coorN = (shotsetting[4].receiver_n[1]:shotsetting[4].receiver_n[2]:shotsetting[4].receiver_n[3]) .-
            shotsetting[4].modelrange[1] / 2.0
    coorE = (shotsetting[4].receiver_e[1]:shotsetting[4].receiver_e[2]:shotsetting[4].receiver_e[3]) .-
            shotsetting[4].modelrange[2] / 2.0
    coorD = (shotsetting[4].receiver_d[1]:shotsetting[4].receiver_d[2]:shotsetting[4].receiver_d[3])
    Ln = length(coorN)
    Le = length(coorE)
    Ld = length(coorD)
    Lt = shotsetting[1].outnpts
    idtable = zeros(Int, 3, Ld * Le * Ln)
    for ix = 1:Le, iy = 1:Ln, iz = 1:Ld
        ir = SEMshell.stationnumber(ix, iy, iz, Le, Ln, Ld)
        idtable[1, ir] = iz
        idtable[2, ir] = ix
        idtable[3, ir] = iy
    end

    Veloc = zeros(Float32, Ld, Le, Ln, 3, Lt)
    DviDxj = zeros(Float32, Ld, Le, Ln, 3, 3, Lt)
    G = zeros(Float32, Ld, Le, Ln, 6, Lt)

    println("Collect data")
    while !all(jobfinish)
        iproc = take!(signchannel)
        lock(lk)
        jid = jobid[iproc-1]
        (ids, velocs, ∇v, g) = fetch(jobbuffer[iproc-1])
        println("Get data for job ", jid, " from worker ", iproc)
        put!(proc_avail, iproc)
        unlock(lk)
        Threads.@threads for irl in eachindex(ids)
            id = ids[irl]
            iz = idtable[1, id]
            iy = idtable[2, id]
            ix = idtable[3, id]
            Veloc[iz, iy, ix, 1, :] .= velocs[2, irl, :]
            Veloc[iz, iy, ix, 2, :] .= velocs[1, irl, :]
            Veloc[iz, iy, ix, 3, :] .= -velocs[3, irl, :]
            DviDxj[iz, iy, ix, 1, 1, :] .= ∇v[2, 2, irl, :]
            DviDxj[iz, iy, ix, 1, 2, :] .= ∇v[1, 2, irl, :]
            DviDxj[iz, iy, ix, 1, 3, :] .= -∇v[3, 2, irl, :]
            DviDxj[iz, iy, ix, 2, 1, :] .= ∇v[2, 1, irl, :]
            DviDxj[iz, iy, ix, 2, 2, :] .= ∇v[1, 1, irl, :]
            DviDxj[iz, iy, ix, 2, 3, :] .= -∇v[3, 1, irl, :]
            DviDxj[iz, iy, ix, 3, 1, :] .= -∇v[2, 3, irl, :]
            DviDxj[iz, iy, ix, 3, 2, :] .= -∇v[1, 3, irl, :]
            DviDxj[iz, iy, ix, 3, 3, :] .= ∇v[3, 3, irl, :]
            G[iz, iy, ix, :, :] .= g[:, irl, :]
        end
        jobfinish[jid] = true
        println("job ", jid, " finish")
    end

    # for i in eachindex(joblist)
    #     (ids, velocs, ∇v, g) = fetch(resultfuture[i])
    #     j = i + input.nproc
    #     if j <= length(joblist)
    #         resultfuture[j] = remotecall(calproc, jobworker[j], joblist[j])
    #     end
    #     Threads.@threads for irl in eachindex(ids)
    #         id = ids[irl]
    #         iz = idtable[1, id]
    #         iy = idtable[2, id]
    #         ix = idtable[3, id]
    #         Veloc[iz, iy, ix, 1, :] .= velocs[2, irl, :]
    #         Veloc[iz, iy, ix, 2, :] .= velocs[1, irl, :]
    #         Veloc[iz, iy, ix, 3, :] .= -velocs[3, irl, :]
    #         DviDxj[iz, iy, ix, 1, 1, :] .= ∇v[2, 2, irl, :]
    #         DviDxj[iz, iy, ix, 1, 2, :] .= ∇v[1, 2, irl, :]
    #         DviDxj[iz, iy, ix, 1, 3, :] .= -∇v[3, 2, irl, :]
    #         DviDxj[iz, iy, ix, 2, 1, :] .= ∇v[2, 1, irl, :]
    #         DviDxj[iz, iy, ix, 2, 2, :] .= ∇v[1, 1, irl, :]
    #         DviDxj[iz, iy, ix, 2, 3, :] .= -∇v[3, 1, irl, :]
    #         DviDxj[iz, iy, ix, 3, 1, :] .= -∇v[2, 3, irl, :]
    #         DviDxj[iz, iy, ix, 3, 2, :] .= -∇v[1, 3, irl, :]
    #         DviDxj[iz, iy, ix, 3, 3, :] .= ∇v[3, 3, irl, :]
    #         G[iz, iy, ix, :, :] .= g[:, irl, :]
    #     end
    # end

    println("Write to disk")
    odir = abspath(input.outputdir)
    mkpath(odir)

    origs = (coorN[1], coorE[1], coorD[1])
    dx = (shotsetting[4].receiver_n[2], shotsetting[4].receiver_e[2], shotsetting[4].receiver_d[2])
    ns = (Ln, Le, Ld)
    cmp = ('x', 'y', 'z')
    Threads.@threads for it in axes(Veloc, 5)
        for i = 1:3
            c = cmp[i]
            writevtkgrid(joinpath(odir, @sprintf("v%c_it%06d.vtk", c, it)), ns, dx, origs, @view(Veloc[:, :, :, i, it]))
            for j = 1:3
                d = cmp[j]
                writevtkgrid(joinpath(odir, @sprintf("H%c%c_it%06d.vtk", c, d, it)), ns, dx, origs,
                             @view(DviDxj[:, :, :, i, j, it]))
            end
        end
        for i = 1:6
            writevtkgrid(joinpath(odir, @sprintf("m%d_it%06d.vtk", i, it)), ns, dx, origs, @view(G[:, :, :, i, it]))
        end
    end
end

main()