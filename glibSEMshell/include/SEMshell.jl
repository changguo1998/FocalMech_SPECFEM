module SEMshell
using Printf, Statistics, LinearAlgebra, FFTW, TOML, SpecialFunctions
import Base.show

const SEMHOME = abspath(@__DIR__, "..", "semexec")
const PRJHOME = abspath(@__DIR__, "../..")
const RERANGELIST = ((1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3))
const HEADER_LEN = 512
const MAX_NETWORK_LEN = 8
const MAX_STATION_LEN = 32
const AMPCORRECTION = 15.17531025
const TAGLEN = MAX_NETWORK_LEN + MAX_STATION_LEN
const DIGITLIST = ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i',
    'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 'A', 'B',
    'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
    'V', 'W', 'X', 'Y', 'Z')
const K = length(DIGITLIST)

include("semmodel.jl")
include("sempartialinterp.jl")

macro must(cond, text="")
    return :(
        if !($(esc(cond)))
            throw($(esc(text)))
        end
    )
end

macro hadbetter(cond, text="")
    return :(
        if !($(esc(cond)))
            @warn($text)
        end
    )
end

macro linearscale(x, xi, xa, yi, ya)
    return :(($(esc(x)) - $(esc(xi))) / ($(esc(xa)) - $(esc(xi))) * ($(esc(ya)) - $(esc(yi))) + $(esc(yi)))
end

# = ============================================================================
# * * * * * * * * * * * *
# *      Settings       *
# * * * * * * * * * * * *
struct ReceiverGroup
    n::Tuple{Float64,Float64,Float64}
    e::Tuple{Float64,Float64,Float64}
    d::Tuple{Float64,Float64,Float64}
    id::Int
end

function ReceiverGroup(n::Tuple{<:Real,<:Real,<:Real},
    e::Tuple{<:Real,<:Real,<:Real},
    d::Tuple{<:Real,<:Real,<:Real},
    id::Integer=1)
    return ReceiverGroup(Float64.(n), Float64.(e), Float64.(d), Int(id))
end

function ReceiverGroup(nstart::Real, nlen::Integer, nstop::Real,
    estart::Real, elen::Integer, estop::Real,
    dstart::Real, dlen::Integer, dstop::Real,
    id::Integer=1)
    return ReceiverGroup((Float64(nstart), Float64((nstop - nstart) / (nlen - 1)), Float64(nstop)),
        (Float64(estart), Float64((estop - estart) / (elen - 1)), Float64(estop)),
        (Float64(dstart), Float64((dstop - nstart) / (dlen - 1)), Float64(dstop)), Int(id))
end

function ReceiverGroup(nstart::Real, Δn::AbstractFloat, nstop::Real,
    estart::Real, Δe::AbstractFloat, estop::Real,
    dstart::Real, Δd::AbstractFloat, dstop::Real,
    id::Integer=1)
    return ReceiverGroup((Float64(nstart), Float64(Δn), Float64(nstop)),
        (Float64(estart), Float64(Δe), Float64(estop)),
        (Float64(dstart), Float64(Δd), Float64(dstop)), Int(id))
end

function show(io::IO, r::ReceiverGroup)
    @printf(io, "%d (%g:%g:%g, %g:%g:%g, %g:%g:%g)", r.id, r.n[1], r.n[2], r.n[3], r.e[1], r.e[2], r.e[3],
        r.d[1], r.d[2], r.d[3])
    return nothing
end

struct ShotSetting
    wkdir::String
    nproc_ξ::Int
    nproc_η::Int
    nex_ξ::Int
    nex_η::Int
    nex_γ::Int
    irregularmesh::Bool
    surfacemesh::Vector{Int}
    step_multiply::Int
    outdt::Float64
    outnpts::Int
    srcloc::Tuple{Float64,Float64,Float64}
    forceDirc::Char
    M0::Float64
    risetime::Float64
    modelrange::Tuple{Float64,Float64,Float64}
    model::String
    pml::Tuple{Float64,Float64,Float64}
    # receiver_n::Tuple{Float64,Float64,Float64}
    # receiver_e::Tuple{Float64,Float64,Float64}
    # receiver_d::Tuple{Float64,Float64,Float64}
    receiver::Vector{ReceiverGroup}
end

function ShotSetting(wkdir::AbstractString, nproc_ξ::Real, nproc_η::Real, nex_ξ::Real, nex_η::Real, nex_γ::Real,
    irregularmesh::Bool, surfacemesh::AbstractVector{<:Integer}, step_multiply::Real, outdt::Real,
    outnpts::Real, srcloc,
    forceDirc::Char, M0::Real, rt::Real, modelrange, model::AbstractString, pml, receiver_n,
    receiver_e, receiver_d, receivertag::Integer=1)
    mr = (Float64(modelrange[1]),
        Float64(modelrange[2]),
        Float64(modelrange[3]))
    _surf = zeros(Int, length(surfacemesh))
    for i in eachindex(surfacemesh)
        _surf[i] = Int(surfacemesh[i])
    end
    _src = (Float64(srcloc[1]), Float64(srcloc[2]), Float64(srcloc[3]))
    pml_ = (Float64(pml[1]), Float64(pml[2]), Float64(pml[3]))
    rn = (Float64(receiver_n[1]), Float64(receiver_n[2]), Float64(receiver_n[3]))
    re = (Float64(receiver_e[1]), Float64(receiver_e[2]), Float64(receiver_e[3]))
    rd = (Float64(receiver_d[1]), Float64(receiver_d[2]), Float64(receiver_d[3]))
    @must forceDirc in ('N', 'E', 'D', 'T') "force direction must be one of 'N', 'E', 'D' or 'T'"
    # return ShotSetting(String(wkdir), Int(nproc_ξ), Int(nproc_η), Int(nex_ξ), Int(nex_η), Int(nex_γ), irregularmesh,
    #     _surf, Int(step_multiply), Float64(outdt), Int(outnpts), Float64(srcdep), forceDirc,
    #     Float64(M0), Float64(rt), mr, String(model), pml_, rn, re, rd)
    return ShotSetting(String(wkdir), Int(nproc_ξ), Int(nproc_η), Int(nex_ξ), Int(nex_η), Int(nex_γ), irregularmesh,
        _surf, Int(step_multiply), Float64(outdt), Int(outnpts), _src, forceDirc,
        Float64(M0), Float64(rt), mr, String(model), pml_, [ReceiverGroup(rn, re, rd, receivertag)])
end

function ShotSetting(wkdir::AbstractString, nproc_ξ::Real, nproc_η::Real, nex_ξ::Real, nex_η::Real, nex_γ::Real,
    irregularmesh::Bool, surfacemesh::AbstractVector{<:Integer}, step_multiply::Real, outdt::Real,
    outnpts::Real, srcloc,
    forceDirc::Char, M0::Real, rt::Real, modelrange, model::AbstractString, pml, rcvs::Vector{ReceiverGroup})
    mr = (Float64(modelrange[1]),
        Float64(modelrange[2]),
        Float64(modelrange[3]))
    _surf = zeros(Int, length(surfacemesh))
    for i in eachindex(surfacemesh)
        _surf[i] = Int(surfacemesh[i])
    end
    _src = (Float64(srcloc[1]), Float64(srcloc[2]), Float64(srcloc[3]))
    pml_ = (Float64(pml[1]), Float64(pml[2]), Float64(pml[3]))
    # rn = (Float64(receiver_n[1]), Float64(receiver_n[2]), Float64(receiver_n[3]))
    # re = (Float64(receiver_e[1]), Float64(receiver_e[2]), Float64(receiver_e[3]))
    # rd = (Float64(receiver_d[1]), Float64(receiver_d[2]), Float64(receiver_d[3]))
    @must forceDirc in ('N', 'E', 'D', 'T') "force direction must be one of 'N', 'E', 'D' or 'T'"
    # return ShotSetting(String(wkdir), Int(nproc_ξ), Int(nproc_η), Int(nex_ξ), Int(nex_η), Int(nex_γ), irregularmesh,
    #     _surf, Int(step_multiply), Float64(outdt), Int(outnpts), Float64(srcdep), forceDirc,
    #     Float64(M0), Float64(rt), mr, String(model), pml_, rn, re, rd)
    return ShotSetting(String(wkdir), Int(nproc_ξ), Int(nproc_η), Int(nex_ξ), Int(nex_η), Int(nex_γ), irregularmesh,
        _surf, Int(step_multiply), Float64(outdt), Int(outnpts), _src, forceDirc,
        Float64(M0), Float64(rt), mr, String(model), pml_, rcvs)
end

function show(io::IO, s::ShotSetting)
    println(io, "wkdir: ", s.wkdir)
    println(io, "model: ", s.model)
    println(io, "time:  ", s.outdt, "s x ", s.outnpts, "; outdt = ", s.outdt/s.step_multiply, "s x ", s.step_multiply)
    println(io, "mesh  (", s.irregularmesh ? ("irregular " * join(s.surfacemesh, ' ')) : "regular", ") :")
    println(io, "  N ", s.modelrange[1], "km + 2 * ", s.pml[1], "km ", s.nex_η, "/", s.nproc_η)
    println(io, "  E ", s.modelrange[2], "km + 2 * ", s.pml[2], "km ", s.nex_ξ, "/", s.nproc_ξ)
    println(io, "  D ", s.modelrange[3], "km + ", s.pml[3], "km ", s.nex_γ)
    println(io, "source:")
    println(io, "  loc  ", s.srcloc)
    println(io, "  dirc ", s.forceDirc)
    println(io, "  M0   ", s.M0)
    println(io, "  rt   ", s.risetime, "s")
    println(io, "receiver:")
    for r in s.receiver
        print(io, "  ")
        show(io, r)
        print(io, '\n')
    end
    return nothing
end

# = ============================================================================
# * * * * * * * * * * * *
# *   Basic functions   *
# * * * * * * * * * * * *

function showsetting(s::ShotSetting)
    t = String["Current Settings:"]
    for varname in keys(s)
        push!(t,
            "\t" * string(varname) * " = " *
            string(getproperty(s, Symbol(varname))))
    end
    push!(t, @sprintf("\tX step: %.2f", s.modelrange[1] / s.nex_ξ))
    push!(t, @sprintf("\tY step: %.2f", s.modelrange[2] / s.nex_η))
    push!(t, @sprintf("\tZ step: %.2f", s.modelrange[3] / s.nex_γ))
    @info join(t, "\n")
    return nothing
end

function checksetting(s::ShotSetting)
    # @must s.receiver_n[1] >= 0.0 "X of some receivers are too small"
    # @must s.receiver_e[1] >= 0.0 "Y of some receivers are too small"
    # @must s.receiver_n[end] <= s.modelrange[1] "X of some receivers are too large"
    # @must s.receiver_e[end] <= s.modelrange[2] "Y of some receivers are too large"
    # @must s.receiver_d[end] <= s.modelrange[3] "Z of some receivers are too large"
    for r in s.receiver
        @must r.n[1] >= 0.0 "In group " * string(r.id) * " X of some receivers are too small"
        @must r.e[1] >= 0.0 "In group " * string(r.id) * " Y of some receivers are too small"
        @must r.n[end] <= s.modelrange[1] "In group " * string(r.id) * " X of some receivers are too large"
        @must r.e[end] <= s.modelrange[2] "In group " * string(r.id) * " Y of some receivers are too large"
        @must r.d[end] <= s.modelrange[3] "In group " * string(r.id) * " Z of some receivers are too large"
    end
    return nothing
end

function multireplace(str::String, pat...)
    t = deepcopy(str)
    for p in pat
        t = replace(t, p)
    end
    return t
end

function cleandir(path::String)
    mkpath(path)
    rm.(readdir(path; join=true), recursive=true, force=true)
    return nothing
end

function lines2mat!(d, s::Vector{T}) where {T<:AbstractString}
    for i in axes(d, 1)
        l = split(s[i]; keepempty=false)
        for j in axes(d, 2)
            d[i, j] = parse(Float32, l[j])
        end
    end
    return nothing
end

_ceil(x::Real, y::Real) = ceil(x / y) * y

_floor(x::Real, y::Real) = floor(x / y) * y

# = ============================================================================
# * * * * * * * * * * * *
# *       Model         *
# * * * * * * * * * * * *

function resamplemodel(s::ShotSetting)
    println("Load model")
    model1 = Model(s.model)
    # * reconstruct material hash
    v0 = Float64.(model1.v0)
    dv = Float64.(model1.dv)
    nv = Int.(model1.np)

    mathash = zeros(Int, nv...)
    for i in eachindex(mathash)
        mathash[i] = i
    end
    (nd, ne, nn) = size(model1.index)
    println("Add PML layer")
    # * append pml layers
    npmln = round(Int, s.pml[1] / s.modelrange[1] * (nn - 1))
    npmle = round(Int, s.pml[2] / s.modelrange[2] * (ne - 1))
    npmld = round(Int, s.pml[3] / s.modelrange[3] * (nd - 1))
    # * vp, vs, rho field
    vp = zeros(nd, ne, nn)
    vs = zeros(nd, ne, nn)
    rho = zeros(nd, ne, nn)
    air = falses(nd, ne, nn)
    Threads.@threads for idx = 1:nn*ne*nd
        mid = model1.index[idx]
        if iszero(mid)
            vp[idx] = 0.34
            vs[idx] = 0.0001
            rho[idx] = 0.00129
            air[idx] = true
        else
            vp[idx] = model1.table[1, mid]
            vs[idx] = model1.table[2, mid]
            rho[idx] = model1.table[3, mid]
            air[idx] = false
        end
    end
    (gridfields, gridair) = appendpml(vp, vs, rho, air, npmln, npmle, npmld)
    # println(minimum(gridfields[1][.!gridair]))
    gridmaterial = field2material(gridfields, v0, dv, mathash, gridair)
    GC.gc()
    println("Construct spec element")
    # * resample mesh
    gridmask = maskgridinterface(gridmaterial)
    semmask = falses(s.nex_γ + 2, s.nex_ξ + 1, s.nex_η + 1)
    semair = falses(s.nex_γ + 1, s.nex_ξ, s.nex_η)
    semair[1, :, :] .= true
    if s.irregularmesh
        semmesh = shiftmesh_onlytopo!(semmask, gridmask)
    else
        semmesh = shiftmesh_split!(semmask, gridmask)
    end
    interpunmaskedcoor!(semmesh, semmask)
    smoothmeshcoor!(@view(semmesh[1, 3:end, :, :]), 1)
    interpunmaskedcoor!(semmesh, semmask)
    if s.irregularmesh
        sstruct = deepcopy(s.surfacemesh)
        fixedlayer = cumsum(sstruct .+ 2)
        factor = fill(Int(2), length(sstruct))
        factor[1] = 1
        cumprod!(factor, factor)
        # for idx in CartesianIndices((s.nex_η + 1, s.nex_ξ + 1))
        el = mean(semmesh[1, 2, :, :])
        # el = semmesh[1, 2, idx]
        # h = (1.0 - el) / (sum(factor .* sstruct) + factor[end] * 2 * (s.nex_γ - sum(sstruct)))
        h = min((s.modelrange[1] + 2 * s.pml[1]) / s.nex_η, (s.modelrange[2] + 2 * s.pml[2]) / s.nex_ξ) /
            (s.modelrange[3] + s.pml[3])
        znorm = zeros(length(sstruct))
        for i in eachindex(znorm)
            if i == 1
                znorm[i] = (sstruct[i]+1)*h + el
            else
                znorm[i] = znorm[i-1] + h * factor[i] * (sstruct[i]+2)
            end
        end
        for i in eachindex(znorm)
            l = fixedlayer[i] + 1
            # l2 = l + 1
            # semmask[l, idx] = true
            semmask[l, :, :] .= true
            # semmask[l2, :, :] .= true
            for idx in CartesianIndices((s.nex_ξ + 1, s.nex_η + 1))
                semmesh[1, l, idx] = znorm[i]
                # semmesh[1, l2, idx] = znorm[i] + 3 * h * factor[i]
            end
        end
        # end
        # semmask[NTHINLAYER+3, :, :] .= true
        # k = h + (NTHINLAYER + 1) * (1.0 - h) / (2 * s.nex_γ - NTHINLAYER - 1)
        # for idx in CartesianIndices((s.nex_η + 1, s.nex_ξ + 1))
        #     semmesh[1, NTHINLAYER+3, idx] = k
        # end
        interpunmaskedcoor!(semmesh, semmask)
    end

    # * interpolate material
    println("Interpolate material")
    # semfields            = interpfield2cell(semmesh, gridfields)
    semfields = interpfield2cell_looponfd(semmesh, gridfields, gridair)
    println("Convert to material")
    semmaterial_oldindex = field2material(semfields, v0, dv, mathash, semair)
    println("Set new index for material")
    oldmats = filter(>(0), unique(semmaterial_oldindex))
    newmattable = model1.table[:, oldmats]
    semmaterial = zeros(eltype(semmaterial_oldindex), size(semmaterial_oldindex))
    Threads.@threads for idx in eachindex(semmaterial_oldindex)
        for i in eachindex(oldmats)
            if oldmats[i] == semmaterial_oldindex[idx]
                semmaterial[idx] = i
                break
            end
        end
    end
    println("finish model")
    return (newmattable, semmaterial, semmesh, model1)
end

_pmlwidth(wref::Real, range::Real, n::Integer) = _floor(wref, range / n) - range / n / 2

function writemeshpar(s::ShotSetting)
    (mtable, tmindex, tmcoor, mdl) = resamplemodel(s); GC.gc()
    mindex = tmindex[2:end, :, :]
    mcoor = tmcoor[:, 2:end, :, :]
    println("Split regions")
    region = cell2region(mindex)
    nlayer = size(mindex, 1)
    h = @. (s.modelrange + (2, 2, 1) * s.pml) / (s.nex_η, s.nex_ξ, s.nex_γ)
    doublinglayer = ""
    if s.irregularmesh
        ss = s.surfacemesh .+ 2
        ss[1] -= 1
        for i in eachindex(s.surfacemesh)
            doublinglayer *= @sprintf("NZ_DOUBLING_%d = %d\n", i, s.nex_γ - sum(ss[1:i]) +1)
        end
    end
    # * Mesh_Par_file
    println("Write Mesh_Par_file")
    pat = read(joinpath(@__DIR__, "Mesh_Par_file.template"), String)
    pat = multireplace(pat,
        "{{LATITUDE_MIN}}" => @sprintf("%g", -1000 * s.pml[1]),
        "{{LATITUDE_MAX}}" => @sprintf("%f", (s.modelrange[1] + s.pml[1]) * 1000),
        "{{LONGITUDE_MIN}}" => @sprintf("%g", -1000 * s.pml[2]),
        "{{LONGITUDE_MAX}}" => @sprintf("%f", (s.modelrange[2] + s.pml[2]) * 1000),
        "{{DEPTH_BLOCK_KM}}" => @sprintf("%.1f", s.modelrange[3] + s.pml[3]),
        "{{NEX_XI}}" => @sprintf("%d", s.nex_ξ),
        "{{NEX_ETA}}" => @sprintf("%d", s.nex_η),
        "{{NPROC_XI}}" => @sprintf("%d", s.nproc_ξ),
        "{{NPROC_ETA}}" => @sprintf("%d", s.nproc_η),
        "{{PML_X_THICKNESS}}" => @sprintf("%g", _pmlwidth(s.pml[2], s.modelrange[2] + 2 * s.pml[2], s.nex_ξ) * 1000),
        "{{PML_Y_THICKNESS}}" => @sprintf("%g", _pmlwidth(s.pml[1], s.modelrange[1] + 2 * s.pml[1], s.nex_η) * 1000),
        "{{PML_Z_THICKNESS}}" => @sprintf("%g", _pmlwidth(s.pml[3], s.modelrange[3] +     s.pml[3], s.nex_γ) * 1000),
        "{{NDOUBLINGS}}" => @sprintf("%d", length(s.surfacemesh)),
        "{{MULTIPLE_LAYER_NUMBER}}" => doublinglayer,
        "{{REGULAR_MESH}}" => s.irregularmesh ? ".false." : ".true.")
    open(abspath(s.wkdir, "DATA", "meshfem3D_files", "Mesh_Par_file"), "w") do io
        println(io, pat)
        println(io, "NMATERIALS = ", size(mtable, 2))
        for i in axes(mtable, 2)
            @printf(io, "%d %g %g %g 999.0 999.0 0 2\n", i, mtable[3, i] * 1000.0,
                mtable[1, i] * 1000.0, mtable[2, i] * 1000.0)
        end
        println(io, "NREGIONS = ", size(region, 2))
        for i in axes(region, 2)
            @printf(io, "%d %d %d %d %d %d %d\n",
                region[3, i], region[4, i],
                region[5, i], region[6, i],
                nlayer - region[2, i] + 1, nlayer - region[1, i] + 1,
                region[7, i])
        end
    end
    # * interfaces.dat
    println("Write interfaces.dat")
    open(abspath(s.wkdir, "DATA", "meshfem3D_files", "interfaces.dat"), "w") do io
        println(io, nlayer)
        for i = 1:nlayer
            @printf(io, ".true. %d %d %g %g %g %g\n",
                size(mcoor, 3), size(mcoor, 4),
                -1000 * s.pml[2], -1000 * s.pml[1],
                1000 * h[2], 1000 * h[1])
            println(io, "interface" * string(i) * ".dat")
        end
        for _ = 1:nlayer
            println(io, "1")
        end
    end

    println("Write interface depth file")
    for l = 1:nlayer
        open(abspath(s.wkdir, "DATA", "meshfem3D_files", "interface" * string(nlayer - l + 1) * ".dat"),
            "w") do io
            for inorth in axes(mcoor, 4), ieast in axes(mcoor, 3)
                # println(io, -1000.0 * mcoor[1, l, ieast, inorth] * (s.modelrange[3] + s.pml[3]))
                println(io,
                    -1000.0 *
                    @linearscale(mcoor[1, l, ieast, inorth], 0.0, 1.0, mdl.zt, (s.modelrange[3] + s.pml[3])))
            end
        end
    end
    return nothing
end

# = ============================================================================
# * * * * * * * * * * * *
# *       Station       *
# * * * * * * * * * * * *

function encodestation(n::Int, L::Int)
    t = n
    s = Char.(zeros(UInt8, L))
    for i = L:-1:1
        (d, r) = divrem(t, K)
        s[i] = DIGITLIST[r+1]
        t = d
    end
    return String(s)
end

function decodestation(s::AbstractString)
    ss = lstrip(s, '0')
    t = 0
    for i in eachindex(ss)
        t *= K
        t += findfirst(DIGITLIST) do x
            x == ss[i]
        end - 1
    end
    return t
end

function stationnumber(ix::Int, iy::Int, iz::Int, Lx::Int, Ly::Int, Lz::Int)
    return (iz - 1) * Lx * Ly + (iy - 1) * Lx + ix
end

@inline function _range_length(a::Real, s::Real, b::Real)
    return floor(Int, (b-a)/s) + 1
end

function stationtable(rcv::Vector{ReceiverGroup})
    nrcv = map(rcv) do r
        # nn = floor(Int, (r.n[3] - r.n[1]) / r.n[2]) + 1
        # ne = floor(Int, (r.e[3] - r.e[1]) / r.e[2]) + 1
        # nd = floor(Int, (r.d[3] - r.d[1]) / r.d[2]) + 1
        nn = _range_length(r.n[1], r.n[2], r.n[3])
        ne = _range_length(r.e[1], r.e[2], r.e[3])
        nd = _range_length(r.d[1], r.d[2], r.d[3])
        nn*ne*nd
    end
    table = zeros(Int, sum(nrcv), 5)
    p = 0
    for r in rcv
        # ys = r.n[1]:r.n[2]:r.n[3]
        # xs = r.e[1]:r.e[2]:r.e[3]
        # zs = r.d[1]:r.d[2]:r.d[3]
        # Lx = length(xs)
        # Ly = length(ys)
        # Lz = length(zs)
        Ly = _range_length(r.n[1], r.n[2], r.n[3])
        Lx = _range_length(r.e[1], r.e[2], r.e[3])
        Lz = _range_length(r.d[1], r.d[2], r.d[3])
        for iy = 1:Ly, ix = 1:Lx, iz = 1:Lz
            p += 1
            table[p, 1] = p
            table[p, 2] = r.id
            table[p, 3] = iz
            table[p, 4] = ix
            table[p, 5] = iy
        end
    end
    return table
end

function writestationposition(s::ShotSetting)
    # ys = s.receiver_n[1]:s.receiver_n[2]:s.receiver_n[3]
    # xs = s.receiver_e[1]:s.receiver_e[2]:s.receiver_e[3]
    # zs = s.receiver_d[1]:s.receiver_d[2]:s.receiver_d[3]
    # L = Int(log(nextpow(K, length(xs) * length(ys) * length(zs))) / log(K))
    # Lx = length(xs)
    # Ly = length(ys)
    # Lz = length(zs)
    # @must (L <= MAX_NETWORK_LEN + MAX_STATION_LEN) "Too many stations to encode"
    stable = stationtable(s.receiver)
    igtable = Dict{Int, Int}()
    for ir in eachindex(s.receiver)
        igtable[s.receiver[ir].id] = ir
    end
    open(abspath(s.wkdir, "DATA", "STATIONS"), "w") do io
        # for ix = 1:Lx, iy = 1:Ly, iz = 1:Lz
        #     str = encodestation(stationnumber(ix, iy, iz, Lx, Ly, Lz), MAX_NETWORK_LEN + MAX_STATION_LEN)
        #     @printf(io, "%s %s %f %f %f %f\n", str[1+MAX_NETWORK_LEN:end], str[1:MAX_NETWORK_LEN],
        #         ys[iy] * 1000, xs[ix] * 1000, -1000 * zs[iz], -1000 * zs[iz])
        # end
        for i = axes(stable, 1)
            str = encodestation(stable[i, 1], MAX_NETWORK_LEN + MAX_STATION_LEN)
            r = s.receiver[igtable[stable[i, 2]]]
            z = r.d[1] + r.d[2]*(stable[i, 3]-1)
            x = r.e[1] + r.e[2]*(stable[i, 4]-1)
            y = r.n[1] + r.n[2]*(stable[i, 5]-1)
            @printf(io, "%s %s %f %f %f %f\n", str[1+MAX_NETWORK_LEN:end], str[1:MAX_NETWORK_LEN],
                y * 1000.0, x * 1000.0, -1000.0 * z, -1000.0 * z)
        end
    end
    return nothing
end

# = ============================================================================
# * * * * * * * * * * * *
# *       Source        *
# * * * * * * * * * * * *

# """
#     stf(t::Real, t0::Real)

# smoothramp function
# """
# function smoothramp(t::Real, t0::Real; eps::Real=1e-3, shiftratio::Real=1.0)
#     p = 2.0 * atanh(1.0 - eps) * (t / t0 - shiftratio)
#     return (tanh(p) + 1.0) / 2.0
# end

# function dsmoothramp(t::Real, t0::Real; eps::Real=1e-3, shiftratio::Real=1.0)
#     p = 2.0 * atanh(1.0 - eps) * (t / t0 - shiftratio)
#     return (1.0 - tanh(p)^2) * atanh(1.0 - eps) / t0
# end

# function ddsmoothramp(t::Real, t0::Real; eps::Real=1e-3, shiftratio::Real=1.0)
#     p = 2.0 * atanh(1.0 - eps) * (t / t0 - shiftratio)
#     return -4.0 * sinh(p) * atanh(1.0 - eps)^2 / cosh(p)^3 / t0^2
# end

# function merf(t::Real, t0::Real)
#     c = π / 0.83255
#     p = c * (t / t0 - 1.5)
#     return (erf(p) + 1.0) / 2.0
# end

function gauss(t::Real, t0::Real)
    tr = sqrt(-log(0.1)) * t0 / π
    p = (t / tr - 8.0)
    return exp(-p^2) / tr / sqrt(pi)
end

# """
# """
# function intricker(t::Real, t0::Real)
#     c = π / 0.83255
#     p = c * (t / t0 - 1.0)
#     return -2 * c^2 * p * exp(-p^2) / (t0^2) / sqrt(pi)
# end

# function ricker(t::Real, t0::Real)
#     c = π / 0.83255
#     p = c * (t / t0 - 1.0)
#     return -2 * c^3 * (1.0 - 2.0 * p^2) * exp(-p^2) / (t0^3) / sqrt(pi)
# end

function writeforcesolution(s::ShotSetting)
    if s.forceDirc == 'N'
        shot = 1
    elseif s.forceDirc == 'E'
        shot = 2
    else
        shot = 3
    end
    open(abspath(s.wkdir, "DATA", "FORCESOLUTION"), "w") do io
        @printf(io, "FORCE 00%d\n", shot)
        println(io, "time shift: 0.0000")
        println(io, "hdur: 0.0")
        @printf(io, "latorUTM: %g\n", s.srcloc[1] * 1000.0)
        @printf(io, "longorUTM: %g\n", s.srcloc[2] * 1000.0)
        @printf(io, "depth: %g\n", s.srcloc[3] * -1000)
        @printf(io, "factor force source: %g\n", s.M0)
        @printf(io, "component dir vect source E: %.1f\n",
            s.forceDirc == 'E' ? 1.0 : 0.0)
        @printf(io, "component dir vect source N: %.1f\n",
            s.forceDirc == 'N' ? 1.0 : 0.0)
        @printf(io, "component dir vect source Z_UP: %.1f\n",
            (s.forceDirc == 'D') || (s.forceDirc == 'T') ? -1.0 : 0.0)
        @printf(io, "%s\n",
            abspath(s.wkdir, "DATA", "EXTERNAL_SOURCE_TIME_FUNCTION"))
    end
    return nothing
end

function writesourcetimefunction(s::ShotSetting)
    t = range(; start=0.0, step=s.outdt / s.step_multiply, length=s.outnpts * s.step_multiply)
    # stf = smoothramp.(t, s.risetime) .* (sqrt(s.M0) * 10^(AMPCORRECTION / 2.0))
    # stf = dsmoothramp.(t, s.risetime) .* (sqrt(s.M0) * 10^(AMPCORRECTION / 2.0))
    stf = gauss.(t, s.risetime) .* s.M0
    # a = stf[1]
    # b = stf[end]
    # stf .-= a
    # stf .*= b / (b - a)
    # @info "source time function npts: $(length(t)), max: $(maximum(stf)), min: $(minimum(stf))"
    open(abspath(s.wkdir, "DATA", "EXTERNAL_SOURCE_TIME_FUNCTION"), "w") do io
        println(io, s.outdt / s.step_multiply)
        for w in stf
            # println(io, w)
            @printf(io, "%12e\n", w)
        end
    end
    return nothing
end

# = ============================================================================
# * * * * * * * * * * * *
# *       Other         *
# * * * * * * * * * * * *

function writemainpar(s::ShotSetting)
    npts = s.step_multiply * s.outnpts
    pat = read(joinpath(@__DIR__, "Par_file.template"), String)
    pat = multireplace(pat,
        "{{DT}}" => @sprintf("%g", s.outdt / s.step_multiply),
        "{{NPTS}}" => @sprintf("%d", npts),
        "{{NPROC}}" => string(s.nproc_ξ * s.nproc_η),
        "{{INNER_MULTIPLY}}" => @sprintf("%d", s.step_multiply),
        "{{INFO_PRINT_STEP}}" => @sprintf("%d", round(Int, npts / 10)),
        "{{NTSTEP_BETWEEN_OUTPUT_SEISMOS}}" => @sprintf("%d", npts + 1),
        "{{USE_BINARY_OUTPUT}}" => ".true.",
        "{{PML_F0}}" => @sprintf("%g", 1.0 / s.risetime))
    open(abspath(s.wkdir, "DATA", "Par_file"), "w") do io
        println(io, pat)
    end
    return nothing
end

function loadshotsetting(path::AbstractString)
    SETTING = TOML.parsefile(path)
    wkdir = abspath(SETTING["wkdir"])
    pml = SETTING["pml"]
    nproc_ξ = SETTING["nproc_xi"]
    nproc_η = SETTING["nproc_eta"]
    # nproc = nproc_ξ * nproc_η
    nex_ξ = nproc_ξ * 8 * SETTING["xi_multiply"]
    nex_η = nproc_η * 8 * SETTING["eta_multiply"]
    modelrange = Tuple(SETTING["modelrange"])
    surfacestructure = SETTING["surface"]
    h0 = min((modelrange[1] + 2 * pml[1]) / nex_η, (modelrange[2] + 2 * pml[2]) / nex_ξ)
    if SETTING["irregular_mesh"]
        @hadbetter nex_ξ >= 2^length(surfacestructure)
        "suggest that xi_multiply not less than $(2^length(surfacestructure))"
        @hadbetter nex_η >= 2^length(surfacestructure)
        "suggest that eta_multiply not less than $(2^length(surfacestructure))"
        factors = 2 .^ (0:length(surfacestructure)-1)
        ss = surfacestructure .+ 2
        ss[1] -= 1
        thicks = factors .* ss .* h0
        sl = sum(thicks)
        nex_γ = ceil(Int, (modelrange[3] + pml[3] - sl) / h0 / (2^length(surfacestructure))) + sum(ss)
    else
        nex_γ = floor(Int, (modelrange[3] + pml[3]) / h0)
    end
    rcv = ReceiverGroup[]
    for t in SETTING["receiver"]
        r = ReceiverGroup(t["n"][1], t["n"][2], t["n"][3],
                        t["e"][1], t["e"][2], t["e"][3],
                        t["d"][1], t["d"][2], t["d"][3], t["id"])
        push!(rcv, r)
    end

    #* === generate model ===
    shotsetting = map(['N', 'E', 'D']) do dirc
        ShotSetting(joinpath(wkdir, "shot" * dirc),
            nproc_ξ,
            nproc_η,
            nex_ξ,
            nex_η,
            nex_γ,
            SETTING["irregular_mesh"],
            surfacestructure,
            SETTING["step_multiply"],
            SETTING["dt"],
            SETTING["npts"],
            SETTING["srcloc"],
            dirc,
            SETTING["amp"],
            SETTING["risetime"],
            SETTING["modelrange"],
            SETTING["modelfile"],
            pml,
            [SETTING["modelrange"][1] / 2.0, 1.0, SETTING["modelrange"][1] / 2.0],
            [SETTING["modelrange"][2] / 2.0, 1.0, SETTING["modelrange"][2] / 2.0],
            [SETTING["modelrange"][3] / 2.0, 1.0, SETTING["modelrange"][3] / 2.0]
        )
    end
    push!(shotsetting,
        ShotSetting(joinpath(wkdir, "shotT"),
            nproc_ξ,
            nproc_η,
            nex_ξ,
            nex_η,
            nex_γ,
            SETTING["irregular_mesh"],
            surfacestructure,
            SETTING["step_multiply"],
            SETTING["dt"],
            1,
            SETTING["srcloc"],
            'T',
            SETTING["amp"],
            SETTING["risetime"],
            SETTING["modelrange"],
            SETTING["modelfile"],
            pml,
            rcv))
    return shotsetting
end

# = ============================================================================
# * * * * * * * * * * * *
# *      Call SEM       *
# * * * * * * * * * * * *

function genSEMdatabase(s::ShotSetting)
    if !isdir(joinpath(s.wkdir, "bin"))
        mkpath(joinpath(s.wkdir, "bin"))
    end
    for exec in ("xmeshfem3D", "xgenerate_databases", "xspecfem3D", "xcombine_vol_data_vtk")
        symlink(abspath(SEMHOME, "bin", exec), abspath(s.wkdir, "bin", exec))
    end
    return nothing
end

function run_xmeshfem3D(s::ShotSetting)
    local nproc = s.nproc_η * s.nproc_ξ
    if nproc == 1
        @info "Shot $(s.forceDirc) Make mesh"
        run(Cmd(`./bin/xmeshfem3D`; dir=s.wkdir))
    else
        @info "Shot $(s.forceDirc) Make mesh on $(nproc) processors"
        run(Cmd(Cmd(["mpirun", "-np", string(nproc), "./bin/xmeshfem3D"]); dir=s.wkdir))
    end
    return nothing
end

function run_xgenerate_databases(s::ShotSetting)
    local nproc = s.nproc_η * s.nproc_ξ
    if nproc == 1
        @info "Shot $(s.forceDirc) Generate database"
        run(Cmd(`./bin/xgenerate_databases`; dir=s.wkdir))
    else
        @info "Shot $(s.forceDirc) Generate database on $nproc processors"
        run(Cmd(Cmd(["mpirun", "-np", string(nproc), "./bin/xgenerate_databases"]); dir=s.wkdir))
    end
    return nothing
end

function run_xspecfem3D(s::ShotSetting)
    local nproc = s.nproc_η * s.nproc_ξ
    if nproc == 1
        @info "Shot $(s.forceDirc) Run solver"
        run(Cmd(`./bin/xspecfem3D`; dir=s.wkdir))
    else
        @info "Shot $(s.forceDirc) Run solver on $nproc processors"
        run(Cmd(Cmd(["mpirun", "-np", string(nproc), "./bin/xspecfem3D"]); dir=s.wkdir))
    end
    return nothing
end

# = ============================================================================
# * * * * * * * * * * * *
# *     Green Lib       *
# * * * * * * * * * * * *
#=
function readshot!(s::ShotSetting, H::AbstractArray, idtable::AbstractMatrix)
    bufFilePath = joinpath(s.wkdir, "OUTPUT_FILES")
    gfile = normpath(bufFilePath, "globinterp.bin")
    pfiles = filter(readdir(bufFilePath; join=true)) do v
        (_, t) = splitdir(v)
        startswith(t, "proc") && endswith(t, "interp.bin")
    end
    (NDIM, NGLOB, NSPEC, NT, NGLL, NREC, spec2glob, rec2spec, l, network, station, hp, νr) = open(readglob, gfile)
    NPROC = length(pfiles)
    for iproc = 1:NPROC
        @info "Shot $(s.forceDirc) Proc $iproc"
        # io = open(pfiles[iproc])
        # (NRECL, lrec➡rec, LGLOB, lglob➡glob, LSPEC, lspec➡spec, ℓ, ξr, ηr, γr, h′,
        #     _, νrl, veloc) = readprocbin_proc(io, NDIM, NGLOB, NT, NGLL)
        # close(io)
        (_, _, LREC, LGLOB, lglob2glob, lrec2rec, ξr, ηr, γr, lveloc) = open(readlocal, pfiles[iproc])
        loop_on_rec = NSPEC > NREC
        # loop_on_rec = true
        ∂v∂xᵢ = zeros(Float_, NDIM, NDIM, LREC)
        if loop_on_rec
            (tmp, tmp1, tmp2, tmp3, glob2lglob) = prepare_∇rec(NDIM, NGLL, LREC, LGLOB, lglob2glob)
        else
            (tmp, tmp1, tmp2, tmp3, LSPEC, glob2lglob, lspec2spec, spec2lspec) = prepare_∇spec(NDIM, NSPEC, NGLL, LREC,
                lrec2rec, LGLOB,
                lglob2glob, NREC,
                rec2spec)
        end
        GC.gc()
        for it = 1:NT
            for idim = 1:NDIM
                if loop_on_rec
                    # ∇andinterp_looprec!(@view(∂v∂xᵢ[idim, :, :]), @view(lveloc[idim, :, it]), NDIM, NSPEC, NGLL, NREC,
                    #                     rec2spec, spec2glob, νr, hp, l, LREC, LGLOB, lrec2rec, lglob2glob, ξr, ηr, γr)
                    ∇andinterp_looprec!(@view(∂v∂xᵢ[idim, :, :]), @view(lveloc[idim, :, it]), NDIM, NSPEC, NGLL, NREC,
                        rec2spec, spec2glob, νr, hp, l, LREC, LGLOB, lrec2rec, ξr, ηr, γr, glob2lglob,
                        tmp, tmp1, tmp2, tmp3)
                else
                    # ∇andinterp_loopspec!(@view(∂v∂xᵢ[idim, :, :]), @view(lveloc[idim, :, it]), NDIM, NSPEC, NGLL, NREC,
                    #                      rec2spec, spec2glob, νr, hp, l, LREC, LGLOB, lrec2rec, lglob2glob, ξr, ηr, γr)
                    ∇andinterp_loopspec!(@view(∂v∂xᵢ[idim, :, :]), @view(lveloc[idim, :, it]), NDIM, NSPEC, NGLL, NREC,
                        rec2spec, spec2glob, νr, hp, l, LREC, LGLOB, lrec2rec, ξr, ηr, γr, glob2lglob,
                        spec2lspec, LSPEC, lspec2spec, tmp, tmp1, tmp2, tmp3)
                end
            end
            Threads.@threads for irl = 1:LREC
                ir = lrec2rec[irl]
                id = network[ir] * station[ir] |> decodestation
                H[it, 1, idtable[1, id], idtable[2, id], idtable[3, id]] = ∂v∂xᵢ[2, 2, irl] # 11
                H[it, 2, idtable[1, id], idtable[2, id], idtable[3, id]] = ∂v∂xᵢ[1, 1, irl] # 22
                H[it, 3, idtable[1, id], idtable[2, id], idtable[3, id]] = ∂v∂xᵢ[3, 3, irl] # 33
                H[it, 4, idtable[1, id], idtable[2, id], idtable[3, id]] = (∂v∂xᵢ[1, 2, irl] + ∂v∂xᵢ[2, 1, irl]) # 12
                H[it, 5, idtable[1, id], idtable[2, id], idtable[3, id]] = -(∂v∂xᵢ[3, 2, irl] + ∂v∂xᵢ[2, 3, irl]) # 13
                H[it, 6, idtable[1, id], idtable[2, id], idtable[3, id]] = -(∂v∂xᵢ[3, 1, irl] + ∂v∂xᵢ[1, 3, irl]) # 23
            end
        end
    end
    return nothing
end

function constructgreenlib(setting::Vector{ShotSetting}, target::AbstractString)
    @must length(setting) == 4 "Shot setting number not correct"
    coorN = (setting[1].receiver_n[1]:setting[1].receiver_n[2]:setting[1].receiver_n[3]) .-
            setting[1].modelrange[1] / 2.0
    coorE = (setting[1].receiver_e[1]:setting[1].receiver_e[2]:setting[1].receiver_e[3]) .-
            setting[1].modelrange[2] / 2.0
    coorD = (setting[1].receiver_d[1]:setting[1].receiver_d[2]:setting[1].receiver_d[3])
    Ln = length(coorN)
    Le = length(coorE)
    Ld = length(coorD)
    Lt = setting[1].outnpts
    idtable = zeros(Int, 3, Ld * Le * Ln)
    for ix = 1:Le, iy = 1:Ln, iz = 1:Ld
        ir = stationnumber(ix, iy, iz, Le, Ln, Ld)
        idtable[1, ir] = iz
        idtable[2, ir] = ix
        idtable[3, ir] = iy
    end

    G = zeros(Float32, Lt, 6, 3, Ld, Le, Ln)
    r = zeros(Int, 4)
    for i = 1:4
        s = setting[i]
        if s.forceDirc == 'N'
            r[1] = i
        elseif s.forceDirc == 'E'
            r[2] = i
        elseif s.forceDirc == 'D'
            r[3] = i
        elseif s.forceDirc == 'T'
            r[4] = i
        end
    end
    @info "Read shots"
    for i = 1:3
        readshot!(setting[r[i]], selectdim(G, 3, i), idtable)
    end
    @info "Write to disk"
    open(joinpath(target, "glib.bin"), "w") do io
        write(io, Float32(setting[1].risetime))
        write(io, Int32(Ln))
        write(io, Int32(Le))
        write(io, Int32(Ld))
        write(io, Int32(Lt))
        write(io, Float32.(coorN))
        write(io, Float32.(coorE))
        write(io, Float32.(coorD))
        write(io, Float32.(range(0.0; step=setting[1].outdt, length=Lt)))
        write(io, G)
    end
    return nothing
end
=#

function readshot!(s::ShotSetting, H::AbstractArray)
    println("read shot: ", s.forceDirc)
    dir1 = joinpath(s.wkdir, "OUTPUT_FILES")
    dir2 = joinpath(String(s.wkdir[1:end-1]*"T"), "OUTPUT_FILES")
    println("config dir: ", dir2, "\nfield dir: ", dir1)
    gfile = normpath(dir2, "globinterp.bin")
    NPROC = length(filter(t->startswith(t, "proc") && endswith(t, "specinfo.bin"), readdir(dir1)))
    (_, _, NSPEC, _, NGLL, NREC, spec2glob, rec2spec, l, network, station, hp, νr) = open(readglob, gfile)
    for iproc = 1:NPROC
        (_, _, LREC, _, lglob2glob, lrec2rec, ξr, ηr, γr) = open(readlocalhead,
            joinpath(dir2, @sprintf("proc%06dlocalinterp.bin", iproc-1)))
        (NDIM, NT, LGLOB, lveloc) = open(readlocalfield,
            joinpath(dir1, @sprintf("proc%06ddispl.bin", iproc-1)))
        id = zeros(Int, LREC)
        Threads.@threads for irl = 1:LREC
            ir = lrec2rec[irl]
            id[irl] = decodestation(network[ir] * station[ir])
        end
        loop_on_rec = NSPEC > NREC
        dvidxj = zeros(Float32, NDIM, NDIM, LREC)
        if loop_on_rec
            (tmp, tmp1, tmp2, tmp3, glob2lglob) = prepare_∇rec(NDIM, NGLL, LREC, LGLOB, lglob2glob)
        else
            (tmp, tmp1, tmp2, tmp3, LSPEC, glob2lglob, lspec2spec, spec2lspec) = prepare_∇spec(NDIM, NSPEC,
                                                                                                        NGLL, LREC,
                                                                                                        lrec2rec,
                                                                                                        LGLOB,
                                                                                                        lglob2glob,
                                                                                                        NREC,
                                                                                                        rec2spec)
        end
        GC.gc()
        for it = 1:NT
            for idim = 1:NDIM
                if loop_on_rec
                    ∇andinterp_looprec!(@view(dvidxj[idim, :, :]), @view(lveloc[idim, :, it]), NDIM,
                                                 NSPEC, NGLL, NREC, rec2spec, spec2glob, νr, hp, l, LREC, LGLOB,
                                                 lrec2rec, ξr, ηr, γr, glob2lglob, tmp, tmp1, tmp2, tmp3)
                else
                    ∇andinterp_loopspec!(@view(dvidxj[idim, :, :]), @view(lveloc[idim, :, it]), NDIM,
                                                  NSPEC, NGLL, NREC, rec2spec, spec2glob, νr, hp, l, LREC, LGLOB,
                                                  lrec2rec, ξr, ηr, γr, glob2lglob, spec2lspec, LSPEC, lspec2spec,
                                                  tmp, tmp1, tmp2, tmp3)
                end
            end # idim
            Threads.@threads for irl = 1:LREC
                ir = lrec2rec[irl]
                sid = decodestation(network[ir] * station[ir])
                H[it, 1, sid] = dvidxj[2, 2, irl] # 11
                H[it, 2, sid] = dvidxj[1, 1, irl] # 22
                H[it, 3, sid] = dvidxj[3, 3, irl] # 33
                H[it, 4, sid] = (dvidxj[1, 2, irl] + dvidxj[2, 1, irl]) # 12
                H[it, 5, sid] = -(dvidxj[3, 2, irl] + dvidxj[2, 3, irl]) # 13
                H[it, 6, sid] = -(dvidxj[3, 1, irl] + dvidxj[1, 3, irl]) # 23
            end # irl
        end # it
    end # iproc
    return nothing
end

function constructgreenlib(setting::Vector{ShotSetting}, target::AbstractString)
    ishott = findfirst(s->s.forceDirc=='T', setting)
    st = setting[ishott]
    println("shot T info:")
    println(st)
    stable = stationtable(st.receiver)
    r = zeros(Int, 4)
    for i = 1:4
        s = setting[i]
        if s.forceDirc == 'N'
            r[1] = i
        elseif s.forceDirc == 'E'
            r[2] = i
        elseif s.forceDirc == 'D'
            r[3] = i
        elseif s.forceDirc == 'T'
            r[4] = i
        end
    end
    Lt = setting[1].outnpts
    G = zeros(Float32, Lt, 6, 3, size(stable, 1))
    @info "Read shots"
    for i = 1:3
        readshot!(setting[r[i]], selectdim(G, 3, i))
    end
    @info "Reshape"
    Hs = Vector{Array{Float32,6}}(undef, length(st.receiver))
    for i = eachindex(st.receiver)
        r = st.receiver[i]
        nn = floor(Int, (r.n[3] - r.n[1]) / r.n[2]) + 1
        ne = floor(Int, (r.e[3] - r.e[1]) / r.e[2]) + 1
        nd = floor(Int, (r.d[3] - r.d[1]) / r.d[2]) + 1
        Hs[i] = zeros(Float32, Lt, 6, 3, nd, ne, nn)
    end
    igtable = Dict{Int,Int}()
    for i in eachindex(st.receiver)
        igtable[st.receiver[i].id] = i
    end
    for i = axes(stable, 1)
        ig = igtable[stable[i, 2]]
        idep = stable[i, 3]
        ieast = stable[i, 4]
        inorth = stable[i, 5]
        Hs[ig][:,:,:,idep,ieast,inorth] .= G[:,:,:,i]
    end
    @info "Write to disk"
    for i = eachindex(Hs)
        r = st.receiver[i]
        coorN = (r.n[1]:r.n[2]:r.n[3]) .- st.srcloc[1]
        coorE = (r.e[1]:r.e[2]:r.e[3]) .- st.srcloc[2]
        coorD = r.d[1]:r.d[2]:r.d[3]
        open(joinpath(target, @sprintf("glib_%d.bin", r.id)), "w") do io
            write(io, Float32(setting[1].risetime))
            write(io, Int32(length(coorN)))
            write(io, Int32(length(coorE)))
            write(io, Int32(length(coorD)))
            write(io, Int32(Lt))
            write(io, Float32.(coorN))
            write(io, Float32.(coorE))
            write(io, Float32.(coorD))
            write(io, Float32.(range(0.0; step=setting[1].outdt, length=Lt)))
            write(io, Hs[i])
        end
    end
    return nothing
end

end
