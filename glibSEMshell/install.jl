#!/usr/bin/env julia

# Compile SEM code

macro runcmd(cmd)
    return :(run(Cmd($cmd, dir=sem_exec_path); wait=true))
end

macro runcmd_noinfo(cmd)
    return :(run(pipeline(Cmd($cmd, dir=sem_exec_path), devnull); wait=true))
end

if !islink("sem-code/external_libs/scotch")
    if isdir("sem-code/external_libs/scotch")
        rm("sem-code/external_libs/scotch"; recursive=true)
    end
    symlink(abspath("sem-code/external_libs/scotch_5.1.12b/"), abspath("sem-code/external_libs/scotch"), dir_target=true)
end

sem_exec_path = joinpath(@__DIR__, "semexec");
if !isdir(sem_exec_path)
    mkpath(sem_exec_path);
else
    rm(sem_exec_path, recursive=true);
    mkpath(sem_exec_path);
end

cmdstr = ["sh", joinpath(@__DIR__, "sem-code", "configure"), "FC=gfortran", "CC=gcc",
    "MPIFC=mpif90", "--with-mpi", "--prefix="*sem_exec_path];
cmd = Cmd(cmdstr);

@runcmd cmd
@runcmd `make xmeshfem3D`
@runcmd `make xgenerate_databases`
@runcmd `make xspecfem3D`
@runcmd `make xcombine_vol_data_vtk`

run(Cmd(`make all`; dir=joinpath(@__DIR__, "slurm/")))
run(Cmd(`julia install.jl`; dir=abspath("NLLOCshell/")))

# @runcmd_noinfo cmd
# @runcmd_noinfo `make xmeshfem3D`
# @runcmd_noinfo `make xgenerate_databases`
# @runcmd_noinfo `make xspecfem3D`
