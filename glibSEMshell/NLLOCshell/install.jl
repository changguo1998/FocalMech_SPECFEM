# * compile NonLinLoc source code

macro runcmd(cmd, dir)
    return :(run(Cmd($cmd, dir=$(esc(dir))); wait=true))
end

nlloc_exec_path = abspath(@__DIR__, "nllocexec")
nlloc_src_path = abspath(@__DIR__, "NonLinLoc/src")
mkpath(nlloc_exec_path)

if islink(abspath(@__DIR__, "NonLinLoc/src/bin"))
    rm(abspath(@__DIR__, "NonLinLoc/src/bin"))
end
symlink(nlloc_exec_path, abspath(@__DIR__, "NonLinLoc/src/bin"))

if isfile(abspath(nlloc_src_path, "CMakeCache.txt"))
    rm(abspath(nlloc_src_path, "CMakeCache.txt"))
end

@runcmd `cmake .` nlloc_src_path
@runcmd `make` nlloc_src_path
