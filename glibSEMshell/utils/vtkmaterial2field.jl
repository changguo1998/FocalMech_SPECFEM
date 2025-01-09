#!/usr/bin/env julia
using ArgumentProcessor

function getvar(vn::AbstractString, lines::AbstractVector, t::Type)
    p = findfirst(contains(vn), lines)
    q = findfirst(==('='), lines[p])
    if t <: AbstractString
        return strip(lines[p][q+1:end])
    else
        return parse(t, strip(lines[p][q+1:end]))
    end
end

addopt!("semshotdir", abbr="S", fmt=" %s", default=" "*pwd())

INP = ArgumentProcessor.parse(ARGS)

semdir = INP.semshotdir

if !isfile(joinpath(semdir, "DATA", "meshfem3D_files", "Mesh_Par_file"))
    error("file not exist: " * joinpath(semdir, "DATA", "meshfem3D_files", "Mesh_Par_file"))
end

semsettings = filter(l->!startswith(l, '#'), readlines(joinpath(semdir, "DATA", "meshfem3D_files", "Mesh_Par_file")))

nmat = getvar("NMATERIALS", semsettings, Int)

mtable = zeros(3, nmat)
for l = semsettings[findfirst(startswith("NMATERIALS"), semsettings) .+ (1:nmat)]
    m = split(l, ' ', keepempty=false)
    mid = parse(Int, m[1])
    rho = parse(Float64, m[2])
    vp = parse(Float64, m[3])
    vs = parse(Float64, m[4])
    mtable[1,mid] = vp/1000
    mtable[2,mid] = vs/1000
    mtable[3,mid] = rho/1000
end

dbdir = joinpath(semdir, "OUTPUT_FILES", "DATABASES_MPI")

meshfiles = filter(contains("_mesh.vtk"), readdir(dbdir))
if isempty(meshfiles)
    error("mesh files not generated")
end

for fn in meshfiles
    lines = readlines(joinpath(dbdir, fn))
    linenum = findfirst(contains("LOOKUP_TABLE"), lines)
    open(joinpath(dbdir, replace(fn, "mesh.vtk"=>"vp.vtk")), "w") do io
        for i = 1:linenum
            println(io, lines[i])
        end
        for i = linenum+1:length(lines)
            if isempty(lines[i])
                continue
            end
            mid = parse(Int, strip(lines[i]))
            println(io, mtable[1,mid])
        end
    end
    open(joinpath(dbdir, replace(fn, "mesh.vtk"=>"vs.vtk")), "w") do io
        for i = 1:linenum
            println(io, lines[i])
        end
        for i = linenum+1:length(lines)
            if isempty(lines[i])
                continue
            end
            mid = parse(Int, strip(lines[i]))
            println(io, mtable[2,mid])
        end
    end
    open(joinpath(dbdir, replace(fn, "mesh.vtk"=>"rho.vtk")), "w") do io
        for i = 1:linenum
            println(io, lines[i])
        end
        for i = linenum+1:length(lines)
            if isempty(lines[i])
                continue
            end
            mid = parse(Int, strip(lines[i]))
            println(io, mtable[3,mid])
        end
    end
end
