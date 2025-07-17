using Pkg

Pkg.activate(".")

Pkg.add("ArgumentProcessor")
Pkg.add("FFTW")

run(Cmd(Cmd(["julia", "install.jl"]); dir=abspath(@__DIR__, "glibSEMshell")))
