#!/usr/bin/env julia
using ArgumentProcessor, TOML

addopt!("position"; abbr="C", fmt=" %s", default=" terminal")
input = ArgumentProcessor.parse(ARGS)

if input.position == "terminal"
    io = stdout
else
    io = open(abspath(input.position), "w")
end

config = Dict("wkdir"          => @__DIR__,
              "modelrange"     => [5.0, 4.0, 3.0],
              "pml"            => [0.5, 0.5, 0.5],
              "srcdep"         => 0.1,
              "dt"             => 0.02,
              "npts"           => 1000,
              "amp"            => 1000.0,
              "risetime"       => 2.0,
              "receiver_n"     => [0.0, 1.0, 5.0],
              "receiver_e"     => [0.0, 1.0, 4.0],
              "receiver_d"     => [0.0, 1.0, 3.0],
              "modelfile"      => abspath(@__DIR__, "model.bin"),
              "nproc_xi"       => 2,
              "nproc_eta"      => 2,
              "ex_factor"      => 4,
              "step_multiply"  => 5,
              "irregular_mesh" => true,
              "surface"        => [3])

TOML.print(io, config, sorted=true)

if input.position != "terminal"
    close(io)
end