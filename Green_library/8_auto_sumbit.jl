SIMULATE = false
printinfo = false

contain_filename(p::String, pat) = readdir(p) |> 
    t->filter(contains(pat), t) |> isempty |> v->!v

stage_wconf(x) = x > 0.0
stage_shotT(x) = x > 1.0
stage_shotA(x) = x >= 2.7
stage_shotD(x) = stage_shotA(x) || (round(UInt8, mod(x*10, 10)) & 0x01) > 0
stage_shotE(x) = stage_shotA(x) || (round(UInt8, mod(x*10, 10)) & 0x02) > 0
stage_shotN(x) = stage_shotA(x) || (round(UInt8, mod(x*10, 10)) & 0x04) > 0
stage_cglib(x) = x > 3.0
stage_cmprs(x) = x > 4.0
stage_fillI(x) = x > 5.0

function test_calculation_stage(path::String)
    fstage = 0.0;
    if isfile(joinpath(path, "writeconfig.flag"))
        fstage = 1.0
    end
    if stage_wconf(fstage) && isfile(joinpath(path, "shotT.flag"))
        fstage = 2.0;
    end
    if stage_shotT(fstage) && isfile(joinpath(path, "shotD.flag"))
        fstage += 0.1;
    end
    if stage_shotT(fstage) && isfile(joinpath(path, "shotE.flag"))
        fstage += 0.2;
    end
    if stage_shotT(fstage) && isfile(joinpath(path, "shotN.flag"))
        fstage += 0.4;
    end
    if round(fstage, digits=1) == 2.7
        fstage = 3.0 # finish all shot
    end
    if stage_shotA(fstage) && isfile(joinpath(path, "calcglib.flag"))
        fstage = 4.0
    end
    if stage_cglib(fstage) && isfile(joinpath(path, "compress.flag"))
        fstage = 5.0
    end
    if stage_cmprs(fstage) && isfile(joinpath(path, "fillinfo.flag"))
        fstage = 6.0
    end

    pstage = fstage;
    if (pstage == 0.0) && any(isdir, joinpath.(path, "shot".*["T", "D", "E", "N"]))
        pstage = 1.0; # write config
    end
    tdir = joinpath(path, "shotT", "OUTPUT_FILES", "DATABASES_MPI")
    if (pstage == 1.0) && isdir(tdir) && contain_filename(tdir, "proc")
        pstage = 2.0; # shotT
    end
    tdir = joinpath(path, "shotD", "OUTPUT_FILES", "DATABASES_MPI")
    if (pstage == 2.0) && stage_shotT(pstage) &&
        isdir(tdir) &&
        contain_filename(tdir, "proc")
        pstage += 0.1; # shotD
    end
    tdir = joinpath(path, "shotE", "OUTPUT_FILES", "DATABASES_MPI")
    if (pstage >= 2.0) && (pstage < 3.0) && stage_shotT(pstage) &&
        isdir(tdir) &&
        contain_filename(tdir, "proc")
        pstage += 0.2; # shotE
    end
    tdir = joinpath(path, "shotN", "OUTPUT_FILES", "DATABASES_MPI")
    if (pstage >= 2.0) && (pstage < 3.0) && stage_shotT(pstage) &&
        isdir(tdir) &&
        contain_filename(tdir, "proc")
        pstage += 0.4; # shotN
    end
    if round(pstage, digits=1) == 2.7
        pstage = 3.0 # finish all shot
    end
    if (pstage == 3.0) && stage_shotA(pstage) && contain_filename(path, r"glib_tmp[0-9]_[0-9]*\.bin")
        pstage = 4.0 # calc glib
    end
    if (pstage == 4.0) && stage_cglib(pstage) && contain_filename(path, r"glib_[0-9]*\.bin")
        pstage = 5.0 # compress glib
    end
    if (pstage == 5.0) && stage_cmprs(pstage) && isdir(joinpath(path, "tlibvar"))
        pstage = 6.0
    end
    return (fstage, pstage)
end

juliapath = "/public/home/guochang/environment/julia/julia-1.8.5/bin/julia"
sbatchpath = "/opt/gridview/slurm/bin/sbatch"

glib_root = abspath("../../dat/glib/")
station_list = readdir() |> x->filter(startswith("batch"), x) |> x->readlines.(x) |>
    x->sort.(x) |> x->vcat(x...)

status_mat = zeros(Int, length(station_list), 10)
for i = eachindex(station_list)
    glibpath = joinpath(glib_root, station_list[i])
    if !isdir(glibpath)
        status_mat[i, 1] = 0
        continue
    elseif isfile(joinpath(glibpath, "slicemodel.flag"))
        status_mat[i, 1] = 2
    else
        status_mat[i, 1] = 1
    end
    (fstage, pstage) = test_calculation_stage(joinpath(glib_root, station_list[i]))
    if stage_wconf(pstage)
        status_mat[i, 2] = 1
    end
    if stage_shotT(pstage)
        status_mat[i, 3] = 1
    end
    if stage_shotD(pstage)
        status_mat[i, 4] = 1
    end
    if stage_shotE(pstage)
        status_mat[i, 5] = 1
    end
    if stage_shotN(pstage)
        status_mat[i, 6] = 1
    end
    if stage_cglib(pstage)
        status_mat[i, 7] = 1
    end
    if stage_cmprs(pstage)
        status_mat[i, 8] = 1
    end
    if stage_fillI(pstage)
        status_mat[i, 9] = 1
    end
    # = = = 
    if stage_wconf(fstage)
        status_mat[i, 2] = 2
    end
    if stage_shotT(fstage)
        status_mat[i, 3] = 2
    end
    if stage_shotD(fstage)
        status_mat[i, 4] = 2
    end
    if stage_shotE(fstage)
        status_mat[i, 5] = 2
    end
    if stage_shotN(fstage)
        status_mat[i, 6] = 2
    end
    if stage_cglib(fstage)
        status_mat[i, 7] = 2
    end
    if stage_cmprs(fstage)
        status_mat[i, 8] = 2
    end
    if stage_fillI(fstage)
        status_mat[i, 9] = 2
    end

    if isfile(joinpath(glibpath, "download.flag"))
        status_mat[i, :] .= 2
    end
end

part_shot = Int[]
part_before_shot = Int[]
part_after_shot = Int[]

global job_pool = 30 - count(isone, status_mat)

global station_stages = zeros(Int, length(station_list))

for i = eachindex(station_list)
    global station_stages
    istage = findfirst(iszero, status_mat[i, :])
    if isnothing(istage)
        station_stages[i] = 10
        continue
    else
        station_stages[i] = istage
    end
    if istage > 6
        push!(part_after_shot, i)
    elseif istage < 4
        push!(part_before_shot, i)
    else
        push!(part_shot, i)
    end
end

part_before_shot = part_before_shot[sortperm(station_stages[part_before_shot], rev=true)]
part_shot = part_shot[sortperm(station_stages[part_shot], rev=true)]
part_after_shot = part_after_shot[sortperm(station_stages[part_after_shot], rev=true)]

if printinfo
    println("job_pool: ", job_pool)
    println("after shot:")
    for i = part_after_shot
        println(station_list[i])
    end
    println("shot:")
    for i = part_shot
        println(station_list[i])
    end
    println("before shot:")
    for i = part_before_shot
        println(station_list[i])
    end
end

for i in part_after_shot
    global job_pool
    glibpath = joinpath(glib_root, station_list[i])
    istage = findfirst(iszero, status_mat[i, :])
    if istage == 7
        if all(==(2), status_mat[i, 4:6]) &&
            (job_pool > 0) &&
            !isfile(joinpath(glibpath, "calcglib.submited"))
            if SIMULATE
                println("submit ", station_list[i], " calcglib")
            else
                touch(joinpath(glibpath, "calcglib.submited"))
                run(Cmd(Cmd([sbatchpath, "calcglib.slurm"]); dir=glibpath))
            end
            job_pool -= 1
        end
    elseif istage == 8
        if (status_mat[i, 7] == 2) &&
            (job_pool > 0) &&
            !isfile(joinpath(glibpath, "compressglib.submited"))
            if SIMULATE
                println("submit ", station_list[i], " compressglib")
            else
                touch(joinpath(glibpath, "compressglib.submited"))
                run(Cmd(Cmd([sbatchpath, "compressglib.slurm"]); dir=glibpath))
            end
            job_pool -= 1;
        end
    elseif istage == 9
        if (status_mat[i, 8] == 2) && (job_pool > 0) && 
            !isfile(joinpath(glibpath, "fillinfo.submited"))
            if SIMULATE
                println("submit ", station_list[i], " fillinfo")
            else
                touch(joinpath(glibpath, "fillinfo.submited"))
                run(Cmd(Cmd([sbatchpath, "fillinfo.slurm"]); dir=glibpath))
            end
            job_pool -= 1;
        end
    end
end

for _i = 1:min(3-length(part_after_shot), length(part_shot))
    global job_pool
    i = part_shot[_i]
    glibpath = joinpath(glib_root, station_list[i])
    istage = findfirst(iszero, status_mat[i, :])
    if status_mat[i, 3] == 2
        run(Cmd(Cmd(["/usr/bin/bash", "copydata.sh"]); dir=glibpath))
        if iszero(status_mat[i, 4]) && (job_pool > 0) && 
            !isfile(joinpath(glibpath, "shotD.submited"))
            if SIMULATE
                println("submit ", station_list[i], " shotD")
            else
                touch(joinpath(glibpath, "shotD.submited"))
                run(Cmd(Cmd([sbatchpath, "submitjob.slurm"]); dir=joinpath(glibpath, "shotD")))
            end
            job_pool -= 1;
        end
        if iszero(status_mat[i, 5]) && (job_pool > 0) && 
            !isfile(joinpath(glibpath, "shotE.submited"))
            if SIMULATE
                println("submit ", station_list[i], " shotE")
            else
                touch(joinpath(glibpath, "shotE.submited"))
                run(Cmd(Cmd([sbatchpath, "submitjob.slurm"]); dir=joinpath(glibpath, "shotE")))
            end
            job_pool -= 1;
        end
        if iszero(status_mat[i, 6]) && (job_pool > 0) && 
            !isfile(joinpath(glibpath, "shotN.submited"))
            if SIMULATE
                println("submit ", station_list[i], " shotN")
            else
                touch(joinpath(glibpath, "shotN.submited"))
                run(Cmd(Cmd([sbatchpath, "submitjob.slurm"]); dir=joinpath(glibpath, "shotN")))
            end
            job_pool -= 1;
        end
    end
end

for i = part_before_shot
    global job_pool
    glibpath = joinpath(glib_root, station_list[i])
    istage = findfirst(iszero, status_mat[i, :])
    if istage == 1
        if job_pool > 0 && !isfile(joinpath(glibpath, "slicemodel.submited"))
            if SIMULATE
                println("submit ", station_list[i], " slicemodel")
            else
                mkpath(glibpath)
                touch(joinpath(glibpath, "slicemodel.submited"))
                run(Cmd(Cmd([sbatchpath, station_list[i]*".slurm"]); dir=joinpath(@__DIR__, "../slurm/")))
            end
            job_pool -= 1;
        end
    elseif istage == 2
        if (status_mat[i, 1] == 2) && (job_pool > 0) && !isfile(joinpath(glibpath, "writeconfig.submited"))
            if SIMULATE
                println("submit ", station_list[i], " writeconfig")
            else
                touch(joinpath(glibpath, "writeconfig.submited"))
                run(Cmd(
                    Cmd([juliapath, "/public/home/guochang/environment/glibSEMshell/slurm/writeslurmscript.jl",
                    "-Cg", "64", "-Tc", "64", "-Tg", "64", "-L", station_list[i]]);
                dir=glibpath))
                run(Cmd(Cmd([sbatchpath, "writeconfig.slurm"]); dir=glibpath))
            end
            job_pool -= 1;
        end
    elseif istage == 3
        if (status_mat[i, 2] == 2) && (job_pool > 0) && !isfile(joinpath(glibpath, "shotT.submited"))
            if SIMULATE
                println("submit ", station_list[i], " shotT")
            else
                touch(joinpath(glibpath, "shotT.submited"))
                run(Cmd(Cmd([sbatchpath, "submitjob.slurm"]); dir=joinpath(glibpath, "shotT")))
            end
            job_pool -= 1;
        end
    end
end
