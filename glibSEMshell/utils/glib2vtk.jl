#!/usr/bin/env julia
using Printf, ArgumentProcessor

include(abspath(@__DIR__, "../include/io.jl"))
# include("../include/io.jl")

function main()
    addopt!("glibfile"; abbr="G", fmt=" %s", required=true, help="Path to green lib file")
    addopt!("vtkdirs"; abbr="V", fmt=" %s", required=true, help="Directory to vtk files")
    addflag!("binary"; abbr="B", help="Use vtk binary format")

    INP = ArgumentProcessor.parse(ARGS)

    (x, y, z, t, H, rt) = glib_readall(INP.glibfile)
    nx = length(x)
    ny = length(y)
    nz = length(z)
    dx = nx == 1 ? 1.0 : x[2] - x[1]
    dy = ny == 1 ? 1.0 : y[2] - y[1]
    dz = nz == 1 ? 1.0 : z[2] - z[1]
    ox = x[1]
    oy = y[1]
    oz = z[1]

    mkpath(INP.vtkdirs)

    c = ('n', 'e', 'd')

    writefunc = INP.binary ? writevtkgrid : printvtkgrid
    for ic = 1:3, im = 1:6
        Threads.@threads for it in eachindex(t)
            vname = @sprintf("%c%d_it%06d.vtk", c[ic], im, it)
            writefunc(abspath(INP.vtkdirs, vname), (nx, ny, nz), (dx, dy, dz), (ox, oy, oz),
                      @view(H[it, im, ic, :, :, :]))
        end
    end
end

main()
