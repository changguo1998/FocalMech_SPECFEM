module NLLOCshell

using Printf

import Base.write

macro must(cond, text="")
    return :(if !($(esc(cond)))
                 throw($(esc(text)))
             end)
end

macro hadbetter(cond, text="")
    return :(if !($(esc(cond)))
                 @warn($text)
             end)
end

const NLLOCEXEC = abspath(@__DIR__, "nllocexec")

# * * * * * * * * * *
# *     Model
# * * * * * * * * * *
struct Model
    v0::Vector{Float32}
    dv::Vector{Float32}
    np::Vector{Int32}
    dx::Float32
    dy::Float32
    dz::Float32
    zt::Float32
    table::Matrix{Float32}
    index::Array{Int32,3}
end

function Model(s::AbstractString)
    io = open(s, "r")
    nf = read(io, Int32)
    v0 = zeros(Float32, nf)
    dv = zeros(Float32, nf)
    np = zeros(Int32, nf)
    read!(io, v0)
    read!(io, dv)
    read!(io, np)
    nm = read(io, Int32)
    table = zeros(Float32, nf, nm)
    read!(io, table)

    nx = read(io, Int32)
    ny = read(io, Int32)
    nz = read(io, Int32)
    dx = read(io, Float32)
    dy = read(io, Float32)
    dz = read(io, Float32)
    zt = read(io, Float32)
    index = zeros(Int32, nz, ny, nx)
    read!(io, index)
    close(io)
    return Model(v0, dv, np, dx, dy, dz, zt, table, index)
end

function write(io::IO, m::Model)
    Base.write(io, Int32(length(m.v0)))
    Base.write(io, Float32.(m.v0))
    Base.write(io, Float32.(m.dv))
    Base.write(io, Int32.(m.np))
    Base.write(io, Int32(size(m.table, 2)))
    Base.write(io, Float32.(m.table))
    Base.write(io, Int32(size(m.index, 3)))
    Base.write(io, Int32(size(m.index, 2)))
    Base.write(io, Int32(size(m.index, 1)))
    Base.write(io, Float32(m.dx))
    Base.write(io, Float32(m.dy))
    Base.write(io, Float32(m.dz))
    Base.write(io, Float32(m.zt))
    Base.write(io, m.index)
    return nothing
end

function resample_nearest(model::Model, h::Real)
    mx = size(model.index, 3)
    my = size(model.index, 2)
    mz = size(model.index, 1)
    nx = floor(Int, (mx - 1)*model.dx/h) + 1
    ny = floor(Int, (my - 1)*model.dy/h) + 1
    nz = floor(Int, (mz - 1)*model.dz/h) + 1
    index = zeros(Int32, nz, ny, nx)
    Threads.@threads for idx in CartesianIndices(index)
        (iz, iy, ix) = idx.I
        jx = max(1, min(mx, round(Int, (ix - 1) * h / model.dx) + 1))
        jy = max(1, min(my, round(Int, (iy - 1) * h / model.dy) + 1))
        jz = max(1, min(mz, round(Int, (iz - 1) * h / model.dz) + 1))
        index[idx] = model.index[jz, jy, jx]
    end
    return Model(model.v0, model.dv, model.np, h, h, h, model.zt, model.table, index)
end

# * * * * * * * * * *
# *     Setting
# * * * * * * * * * *

abstract type NLLocSetting <: Any end

struct TransformType <: NLLocSetting
    name::String
    refEllipsoid::String
    para::Vector{Float64}
end

"""
TransformType(; name="NONE", refEllipsoid="", para::AbstractVector{<:Real})
"""
function TransformType(; name::AbstractString="NONE", refEllipsoid::AbstractString="",
                       para::AbstractVector{<:Real}=Float64[])
    @must uppercase(String(name)) in
          ("GLOBAL", "SIMPLE", "NONE", "SDC", "LAMBERT", "TRANS_MERC", "AZIMUTHAL_EQDIST")
    """transformation type must be one of: "GLOBAL", "SIMPLE", "NONE", "SDC", "LAMBERT", "TRANS_MERC" or "AZIMUTHAL_EQDIST" """
    if name in ("LAMBERT", "TRANS_MERC")
        @must String(refEllipsoid) in
              ("WGS-84", "GRS-80", "WGS-72", "Australian", "Krasovsky", "International", "Hayford-1909", "Clarke-1880",
               "Clarke-1866", "Airy", "Bessel", "Hayford-1830", "Sphere") """Reference ellipsoid must be one of:
               "WGS-84", "GRS-80", "WGS-72", "Australian", "Krasovsky", "International", "Hayford-1909", "Clarke-1880",
               "Clarke-1866", "Airy", "Bessel", "Hayford-1830", "Sphere" """
    end
    return TransformType(uppercase(String(name)), String(refEllipsoid), Float64.(para))
end

function printsetting(io::IO, s::TransformType)
    print(io, "TRANS ", s.name)
    if s.name in ("SIMPLE", "SDC")
        println(io, join(string.(s.para), ' '))
    elseif s.name in ("LAMBERT", "TRANS_MERC")
        println(io, s.name, ' ', join(string.(s.para), ' '))
    else
        println(io, "")
    end
end

struct NLLocSettingGlobal <: NLLocSetting
    control::NamedTuple{(:messagelevel, :randseed),Tuple{Int,Int}}
    trans::TransformType
end

"""
NLLocSettingGlobal(; messagelevel=1, randseed=rand(Int), trans=TransformType())
"""
function NLLocSettingGlobal(; messagelevel::Integer=1, randseed::Integer=mod(rand(Int), 100000),
                            trans::TransformType=TransformType())
    @must messagelevel >= -1 "MessageLevel >= -1"
    return NLLocSettingGlobal((messagelevel=Int(messagelevel), randseed=Int(randseed)), trans)
end

function printsetting(io::IO, s::NLLocSettingGlobal)
    println(io, join(["CONTROL", string(s.control.messagelevel), string(s.control.randseed)], ' '))
    printsetting(io, s.trans)
end

struct NLLocSettingVel2Grid3D <: NLLocSetting
    inputpath::String
    inputfiletype::String
    inputpara::NamedTuple{(:nx, :ny, :nz, :ox, :oy, :oz, :dx, :dy, :dz),
                          Tuple{Int,Int,Int,Float64,Float64,Float64,Float64,Float64,Float64}}
    output::String
    wavetype::Vector{String}
    grid::NamedTuple{(:nx, :ny, :nz, :ox, :oy, :oz, :dx, :dy, :dz, :type),
                     Tuple{Int,Int,Int,Float64,Float64,Float64,Float64,Float64,Float64,String}}
    clip::Tuple{Float64,Float64}
end

"""
NLLocSettingVel2Grid3D(; inputpath="", inputfiletype="",
nxi=2, nyi=2, nzi=2, oxi=0.0, oyi=0.0, ozi=0.0, dxi=1.0, dyi=1.0, dzi=1.0,
output="", wavetypes=String[],
nx=2, ny=2, nz=2, ox=0.0, oy=0.0, oz=0.0, dx=1.0, dy=1.0, dz=1.0,
type="", vmax=0.0, vmin=0.0)
"""
function NLLocSettingVel2Grid3D(; inputpath::AbstractString="", inputfiletype::AbstractString="",
                                nxi::Integer=2, nyi::Integer=2, nzi::Integer=2,
                                oxi::Real=0.0, oyi::Real=0.0, ozi::Real=0.0,
                                dxi::Real=1.0, dyi::Real=1.0, dzi::Real=1.0,
                                output::AbstractString="", wavetypes::Vector{<:AbstractString}=String[],
                                nx::Integer=2, ny::Integer=2, nz::Integer=2,
                                ox::Real=0.0, oy::Real=0.0, oz::Real=0.0,
                                dx::Real=1.0, dy::Real=1.0, dz::Real=1.0,
                                type::AbstractString="", vmax::Real=0.0, vmin::Real=0.0)
    @must vmin < vmax "vmin $(vmin) must be smaller than vmax $(vmax)"
    @must uppercase(String(inputfiletype)) in ("SIMUL2K", "FDTOMO") "input file type must be one of \"SIMUL2K\" or \"FDTOMO\""
    for wt in wavetypes
        @must uppercase(String(wt)) in ("P", "S") "Wave type must be P or S"
    end
    @must uppercase(String(type)) in
          ("VELOCITY", "VELOCITY_METERS", "SLOWNESS", "VEL2", "SLOW2", "SLOW_2_METERS", "SLOW_LEN")
    inppara = (nx=Int(nxi), ny=Int(nyi), nz=Int(nzi), ox=Float64(oxi), oy=Float64(oyi), oz=Float64(ozi),
               dx=Float64(dxi), dy=Float64(dyi), dz=Float64(dzi))
    return NLLocSettingVel2Grid3D(String(inputpath), uppercase(String(inputfiletype)), inppara, String(output),
                                  uppercase.(String.(wavetypes)),
                                  (nx=Int(nx), ny=Int(ny), nz=Int(nz), ox=Float64(ox), oy=Float64(oy),
                                   oz=Float64(oz), dx=Float64(dx), dy=Float64(dy), dz=Float64(dz), type=String(type)),
                                  (Float64(vmin), Float64(vmax)))
end

function printsetting(io::IO, s::NLLocSettingVel2Grid3D)
    println(io, join(["VGOUT", s.output], ' '))
    for wt in s.wavetype
        println(io, join(["VGTYPE", wt], ' '))
    end
    println(io,
            join(["VGGRID",
                  @sprintf("%d %d %d %g %g %g %g %g %g", s.grid.nx, s.grid.ny, s.grid.nz,
                           s.grid.ox, s.grid.oy, s.grid.oz, s.grid.dx, s.grid.dy, s.grid.dz),
                  s.grid.type],
                 ' '))
    print(io, join(["VGINP", s.inputpath, s.inputfiletype], ' '))
    if s.inputfiletype == "FDTOMO"
        @printf(io, " %g %g %g %d %d %d %g %g %g\n", s.inputpara.ox, s.inputpara.oy, s.inputpara.oz,
                s.inputpara.nx, s.inputpara.ny, s.inputpara.nz, s.inputpara.dx, s.inputpara.dy, s.inputpara.dz)
    else
        println(io, "")
    end
    println(io, join(["VGCLIP", string(s.clip[1]), string(s.clip[2])], ' '))
end

struct NLLocSettingSourceLocation <: NLLocSetting
    label::String
    coortype::String
    zsrce::Float64
    elev::Float64
    paraf::Vector{Float64}
    paras::Vector{String}
end

"""
NLLocSettingSourceLocation(; label="", coordinate="", paras=String[], paraf=Float64[], z=0.0, elev=0.0)
"""
function NLLocSettingSourceLocation(; label::AbstractString="", coordinate::AbstractString="",
                                    paras::Vector{<:AbstractString}=String[], paraf::Vector{<:Real}=Float64[],
                                    z::Real=0.0, elev::Real=0.0)
    c = uppercase(String(coordinate))
    @must c in ("XYZ", "LATLON", "LATLONDM", "LATLONDS") "Coordinate type must be one of
    \"XYZ\", \"LATLON\", \"LATLONDM\", \"LATLONDS\""
    if c in ("XYZ", "LATLON")
        @must length(paraf) == 2
    elseif c == "LATLONDM"
        @must length(paraf) == 4
    elseif c == "LATLONDS"
        @must length(paraf) == 6
    end
    if c in ("LATLONDM", "LATLONDS")
        @must length(paras) == 2
        @must uppercase(String(paras[1])) in ("N", "S") "latDir must be N or S"
        @must uppercase(String(paras[2])) in ("W", "E") "lonDir must be W or E"
    end
    return NLLocSettingSourceLocation(String(label), uppercase(String(coordinate)), Float64(z), Float64(elev),
                                      Float64.(paraf), uppercase.(String.(paras)))
end

function printsetting(io::IO, s::NLLocSettingSourceLocation)
    print(io, "GTSRCE ", s.label, " ", s.coortype, " ")
    if s.coortype in ("XYZ", "LATLON")
        print(io, join(string.(s.paraf), ' '), " ")
    elseif c == "LATLONDM"
        print(io, join(string.(s.paraf[1:2]), ' '))
        print(io, s.paras[1], " ")
        print(io, join(string.(s.paraf[3:4]), ' '))
        print(io, s.paras[2], " ")
    elseif c == "LATLONDS"
        print(io, join(string.(s.paraf[1:3]), ' '))
        print(io, s.paras[1], " ")
        print(io, join(string.(s.paraf[4:6]), ' '))
        print(io, s.paras[2], " ")
    end
    println(io, s.zsrce, " ", s.elev)
end

struct NLLocSettingGrid2Time <: NLLocSetting
    gridfile::String
    traveltimefile::String
    wavetype::Vector{String}
    swapbyte::Bool
    gridmode::String
    anglemode::String
    plfd::Tuple{Float64,Int}
    src::Vector{NLLocSettingSourceLocation}
end

"""
NLLocSettingGrid2Time(; gridfile="", ttfile="", wavetype=String[], swapbyte=false, gridmode="",
anglemode="", fdhsinit=0.001, message=1, source=NLLocSettingSourceLocation[])
"""
function NLLocSettingGrid2Time(; gridfile::AbstractString="", ttfile::AbstractString="",
                               wavetypes::Vector{<:AbstractString}=String[], swapbyte::Bool=false,
                               gridmode::AbstractString="", anglemode::AbstractString="",
                               fdhsinit::Real=0.001, message::Integer=1,
                               source::Vector{NLLocSettingSourceLocation}=NLLocSettingSourceLocation[])
    for wt in wavetypes
        @must uppercase(String(wt)) in ("P", "S") "Wave type must be P or S"
    end
    @must uppercase(String(gridmode)) in ("GRID3D", "GRID2D") "Grid type must be GRID3D or GRID2D"
    @must uppercase(String(anglemode)) in ("ANGLES_YES", "ANGLES_NO") "Angle mode must be ANGLES_YES or ANGLES_NO"
    return NLLocSettingGrid2Time(String(gridfile), String(ttfile), uppercase.(String.(wavetypes)), swapbyte,
                                 uppercase(String(gridmode)), uppercase(String(anglemode)),
                                 (Float64(fdhsinit), Int(message)), source)
end

function printsetting(io::IO, s::NLLocSettingGrid2Time)
    for wt in s.wavetype
        println(io, join(["GTFILES", s.gridfile, s.traveltimefile, wt, s.swapbyte ? "1" : "0"], ' '))
    end
    println(io, join(["GTMODE", s.gridmode, s.anglemode], ' '))
    for sr in s.src
        printsetting(io, sr)
    end
    println(io, join(["GT_PLFD", string(s.plfd[1]), string(s.plfd[2])], ' '))
end

function printsetting(io::IO, s::Vector{NLLocSetting})
    for si in s
        printsetting(io, si)
    end
end

function writemodel(io::IO, m::Model, material::Int)
    for iz in axes(m.index, 1), ix in axes(m.index, 3), iy in axes(m.index, 2)
        im = m.index[iz, iy, ix]
        if iszero(im)
            v = 0.01
        else
            v = m.table[material, im]
        end
        write(io, Float32(v))
    end
end

function _swapdim23(x::Array{<:Real, 3})
    nz = size(x, 1)
    ny = size(x, 2)
    nx = size(x, 3)
    y = zeros(eltype(x), nz, nx, ny)
    Threads.@threads for idx in CartesianIndices(y)
        (i, j, k) = idx.I
        y[idx] = x[i, k, j]
    end
    return y
end

function readgrid(root::AbstractString)
    buffer = Pair{Symbol,Any}[]
    formats = readlines(root * ".hdr")
    l1 = split(formats[1]; keepempty=false)
    nx = parse(Int, l1[1])
    ny = parse(Int, l1[2])
    nz = parse(Int, l1[3])
    ox = parse(Float64, l1[4])
    oy = parse(Float64, l1[5])
    oz = parse(Float64, l1[6])
    push!(buffer, :orig => (ox, oy, oz))
    dx = parse(Float64, l1[7])
    dy = parse(Float64, l1[8])
    dz = parse(Float64, l1[9])
    push!(buffer, :delta => (dx, dy, dz))
    type = String(l1[10])
    push!(buffer, :type => type)
    vtype = uppercase(String(l1[11]))
    l2 = split(formats[2]; keepempty=false)
    if type in ("TIME", "TIME2D", "ANGLE", "ANGLE2D")
        elabel = String(l2[1])
        srcx = parse(Float64, l2[2])
        srcy = parse(Float64, l2[3])
        srcz = parse(Float64, l2[4])
        push!(buffer, :src => (elabel, srcx, srcy, srcz))
        l3 = split(formats[3]; keepempty=false)
        trans = String(l3[2])
    else
        trans = String(l2[2])
    end
    push!(buffer, :trans => trans)
    if vtype == "CHAR"
        T = Char
    elseif vtype == "INT"
        T = Int32
    elseif vtype == "SHORT"
        T = Int16
    elseif vtype == "LONG"
        T = Int64
    elseif vtype == "FLOAT"
        T = Float32
    elseif vtype == "DOUBLE"
        T = Float64
    else
        error("invalid data type")
    end
    dat = zeros(T, nz, ny, nx)
    open(root * ".buf", "r") do io
        read!(io, dat)
    end
    push!(buffer, :data => _swapdim23(dat))
    return NamedTuple(buffer)
end

function _read_vector(io::IO, T::Type, n::Integer)
    t = zeros(T, n)
    read!(io, t)
    return t
end

"""
ttlib_readhead(io::IO) -> (nx, ny, nz, dx, dy, dz, ox, oy, oz)
"""
function ttlib_readhead(io::IO)
    seekstart(io)
    (nx, ny, nz) = _read_vector(io, Int32, 3)
    (dx, dy, dz, ox, oy, oz) = _read_vector(io, Float32, 6)
    return (nx, ny, nz, dx, dy, dz, ox, oy, oz)
end

"""
ttlib_readall(io::IO) -> (dx, dy, dz, ox, oy, oz, tt)
"""
function ttlib_readall(io::IO)
    (nx, ny, nz, dx, dy, dz, ox, oy, oz) = ttlib_readhead(io)
    tt = zeros(Float32, 2, nz, ny, nx)
    read!(io, tt)
    return (dx, dy, dz, ox, oy, oz, tt)
end

"""
ttlib_readlocation(io::IO, x::Real, y::Real, z::Real) -> (tp, ts)
"""
function ttlib_readlocation(io::IO, x::Real, y::Real, z::Real)
    (nx, ny, nz, dx, dy, dz, ox, oy, oz) = ttlib_readhead(io)
    ix = floor(Int, (x - ox) / dx) + 1
    iy = floor(Int, (y - oy) / dy) + 1
    iz = floor(Int, (z - oz) / dz) + 1
    if (ix < 1) || (iy < 1) || (iz < 1) || (ix > nx) || (iy > ny) || (iz > nz)
        error("Location out of range")
    end
    if ix == nx
        ix -= 1
    end
    if iy == ny
        iy -= 1
    end
    if iz == nz
        iz -= 1
    end
    h = (x - ox) / dx - ix + 1.0
    k = (y - oy) / dy - iy + 1.0
    l = (z - oz) / dz - iz + 1.0
    seek(io, ((ix - 1) * 2 * nz * ny + (iy - 1) * 2 * nz + (iz - 1) * 2 + 1 - 1 + 9)*4)
    (p000, s000) = _read_vector(io, Float32, 2)
    seek(io, (ix * 2 * nz * ny + (iy - 1) * 2 * nz + (iz - 1) * 2 + 1 - 1 + 9)*4)
    (p100, s100) = _read_vector(io, Float32, 2)
    seek(io, ((ix - 1) * 2 * nz * ny + iy * 2 * nz + (iz - 1) * 2 + 1 - 1 + 9)*4)
    (p010, s010) = _read_vector(io, Float32, 2)
    seek(io, (ix * 2 * nz * ny + iy * 2 * nz + (iz - 1) * 2 + 1 - 1 + 9)*4)
    (p110, s110) = _read_vector(io, Float32, 2)
    seek(io, ((ix - 1) * 2 * nz * ny + (iy - 1) * 2 * nz + iz * 2 + 1 - 1 + 9)*4)
    (p001, s001) = _read_vector(io, Float32, 2)
    seek(io, (ix * 2 * nz * ny + (iy - 1) * 2 * nz + iz * 2 + 1 - 1 + 9)*4)
    (p101, s101) = _read_vector(io, Float32, 2)
    seek(io, ((ix - 1) * 2 * nz * ny + iy * 2 * nz + iz * 2 + 1 - 1 + 9)*4)
    (p011, s011) = _read_vector(io, Float32, 2)
    seek(io, (ix * 2 * nz * ny + iy * 2 * nz + iz * 2 + 1 - 1 + 9)*4)
    (p111, s111) = _read_vector(io, Float32, 2)
    tp = p000 * (1.0 - h) * (1.0 - k) * (1.0 - l) +
        p100 * h * (1.0 - k) * (1.0 - l) +
        p010 * (1.0 - h) * k * (1.0 - l) +
        p110 * h * k * (1.0 - l) +
        p001 * (1.0 - h) * (1.0 - k) * l +
        p101 * h * (1.0 - k) * l +
        p011 * (1.0 - h) * k * l +
        p111 * h * k * l
    ts = s000 * (1.0 - h) * (1.0 - k) * (1.0 - l) +
        s100 * h * (1.0 - k) * (1.0 - l) +
        s010 * (1.0 - h) * k * (1.0 - l) +
        s110 * h * k * (1.0 - l) +
        s001 * (1.0 - h) * (1.0 - k) * l +
        s101 * h * (1.0 - k) * l +
        s011 * (1.0 - h) * k * l +
        s111 * h * k * l
    return (tp, ts)
end

end
