using Mmap, Printf

include("semmodel.jl")

TINT = Union{Signed,Unsigned}
#=
# * * * * * * * * * * * *
# *   FD3D shotframe    *
# * * * * * * * * * * * *

"""
    fd3d_readbin!(fname, H, nx, ny, nz, nt)
"""
function fd3d_readbin!(fname::AbstractString, H::AbstractArray, nx::Integer, ny::Integer, nz::Integer, nt::Integer)
    io = open(fname)
    nr = read(io, Float32) |> Int
    if nr != nx * ny * nz
        error("Number of receiver not match")
    end
    coor = zeros(Float32, 3, nr)
    read!(io, coor)
    w = zeros(Float32, nr, nt)
    read!(io, w)
    close(io)
    idx = floor.(Int,
                 (coor .- repeat(collect(minimum(coor; dims=2)), 1, size(coor, 2))) ./
                 repeat(collect(maximum(coor; dims=2) .- minimum(coor; dims=2)) ./ [nx - 1, ny - 1, nz - 1], 1,
                        size(coor, 2))) .+ 1
    for ir = 1:nr
        H[:, idx[[3, 2, 1], ir]...] .= w[ir, :]
    end
    return nothing
end

"""
    fd3d_readbin(fname, nt) -> H
"""
function fd3d_readbin(fname::AbstractString, nt::TINT)
    io = open(fname, "r")
    nr = read(io, Float32) |> Int
    coor = zeros(Float32, 3, nr)
    read!(io, coor)
    w = zeros(Float32, nr, nt)
    read!(io, w)
    close(io)
    xr = coor[1, :]
    yr = coor[2, :]
    zr = coor[3, :]
    x = unique(xr) |> sort
    nx = length(x)
    dx = x[2] - x[1]
    mx = x[1]
    y = unique(yr) |> sort
    ny = length(y)
    dy = y[2] - y[1]
    my = y[1]
    z = unique(zr) |> sort
    nz = length(z)
    dz = z[2] - z[1]
    mz = z[1]

    H = zeros(Float32, nt, nz, ny, nx)
    for ir = 1:nr
        ix = (xr[ir] - mx) / dx + 1 |> Int
        iy = (yr[ir] - my) / dy + 1 |> Int
        iz = (zr[ir] - mz) / dz + 1 |> Int
        H[:, iz, iy, ix] .= w[ir, :]
    end
    return H
end
=#
# * * * * * * * * * * * *
# *   SEM  shotframe    *
# * * * * * * * * * * * *

const HEADER_LEN = 512
const MAX_NETWORK_LEN = 8
const MAX_STATION_LEN = 32
const TAGLEN = MAX_NETWORK_LEN + MAX_STATION_LEN
const DIGITLIST = ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
    'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 'A', 'B', 'C', 'D',
    'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
    'Y', 'Z')
const K = length(DIGITLIST)

function sem_decodestation(s::AbstractString)
    ss = lstrip(s, '0')
    t = 0
    for i in eachindex(ss)
        t *= K
        t += findfirst(==(ss[i]), DIGITLIST) - 1
    end
    return t
end

function sem_id2xyz(n::TINT, nx::TINT, ny::TINT, nz::TINT)
    ix = mod(n - 1, nx) + 1 |> Int
    iy = mod((n - ix) / nx, ny) + 1 |> Int
    iz = (n - ix - (iy - 1) * nx) / (nx * ny) + 1 |> Int
    return (ix, iy, iz)
end

"""
    sem_readbin!(fname::AbstractString, H::AbstractArray)
"""
function sem_readbin!(fname::AbstractString, H::AbstractArray)
    ts = zeros(UInt8, HEADER_LEN)
    nn = size(H, 4)
    ne = size(H, 3)
    nd = size(H, 2)
    open(fname, "r") do io
        for _ = 1:nn*ne*nd
            read!(io, ts)
            # read!(io, td)
            l = filter(Char.(ts)) do c
                isdigit(c) || isletter(c) || c == '.'
            end
            is = l[[1:MAX_NETWORK_LEN; MAX_NETWORK_LEN+2:TAGLEN+1]] |> String |> sem_decodestation
            (ie_, in_, id_) = sem_id2xyz(is, ne, nn, nd)
            # H[:, id, ie, in] .= td
            read!(io, @view(H[:, id_, ie_, in_]))
        end
        if !eof(io)
            error("Extra data exist")
        end
    end
    return nothing
end

"""
    sem_readbin(fname, nn, ne, nd, nt) -> H
"""
function sem_readbin(fname::AbstractString, nn::TINT, ne::TINT, nd::TINT, nt::TINT)
    H = zeros(Float32, nt, nd, ne, nn)
    sem_readbin!(fname, H)
    return H
end

"""
sem_writebin(fpath, w::AbstractMatrix{Float32}, network, station)
"""
function sem_writebin(fpath::AbstractString, w::AbstractMatrix{Float32},
    network::Vector{String}, station::Vector{String})
    NT = size(w, 1)
    io = open(fpath, "w")
    buf = zeros(UInt8, HEADER_LEN)
    buf[MAX_NETWORK_LEN+1] = UInt8('.')
    buf[TAGLEN+2:TAGLEN+10] = UInt8.(collect(".BHX.semv"))
    for ir = 1:length(station)
        buf[1:MAX_NETWORK_LEN] .= UInt8.(collect(network[ir]))
        buf[MAX_NETWORK_LEN+2:TAGLEN+1] .= UInt8.(collect(station[ir]))
        write(io, buf)
        for it = 1:NT
            write(io, w[it, ir])
        end
    end
    close(io)
    return nothing
end

# * * * * * * * * * * * *
# *         VTK         *
# * * * * * * * * * * * *

function scanvtkgrid(fpath::AbstractString)
    lines = readlines(fpath)
    encodefmt = strip(lines[3]) |> String
    t = split(lines[4])
    datastructure = t[2] |> String
    t = split(lines[5])
    Nindims = parse.(Int, t[2:4])
    t = split(lines[6])
    origin = parse.(Float64, t[2:4])
    t = split(lines[7])
    dx = parse.(Float64, t[2:4])
    t = split(lines[8])
    ndata = parse(Int, t[2])
    t = split(lines[9])
    vname = t[2] |> String
    vtype = t[3] |> String
    fdata = parse.(Float64, lines[11:ndata+10])
    if !((encodefmt == "ASCII") && (datastructure == "STRUCTURED_POINTS") && (vtype == "float"))
        error("file format not correct")
    end
    x = range(; start=origin[1], step=dx[1], length=Nindims[1])
    y = range(; start=origin[2], step=dx[2], length=Nindims[2])
    z = range(; start=origin[3], step=dx[3], length=Nindims[3])
    v = reshape(fdata, Tuple(Nindims))
    return (x, y, z, v)
end

function printvtkgrid(fpath::AbstractString, n::Tuple{TINT,TINT,TINT}, d::Tuple{Real,Real,Real},
    o::Tuple{Real,Real,Real}, H::AbstractArray{T,3}) where {T<:Real}
    open(fpath, "w") do io
        println(io, "# vtk DataFile Version 3.1")
        println(io, "material model VTK file")
        println(io, "ASCII\nDATASET STRUCTURED_POINTS")
        @printf(io, "DIMENSIONS %d %d %d\n", n...)
        @printf(io, "ORIGIN %f %f %f\n", o...)
        @printf(io, "SPACING %f %f %f\n", d...)
        @printf(io, "POINT_DATA %d\n", prod(n))
        println(io, "SCALARS Value float")
        println(io, "LOOKUP_TABLE default")
        for iz = 1:n[3], iy = 1:n[2], ix = 1:n[1]
            println(io, H[iz, iy, ix])
        end
    end
end

function writevtkgrid(fpath::AbstractString, n::Tuple{TINT,TINT,TINT}, d::Tuple{Real,Real,Real},
    o::Tuple{Real,Real,Real}, H::AbstractArray{T,3}) where {T<:Real}
    open(fpath, "w") do io
        println(io, "# vtk DataFile Version 3.1")
        println(io, "material model VTK file")
        println(io, "BINARY\nDATASET STRUCTURED_POINTS")
        @printf(io, "DIMENSIONS %d %d %d\n", n...)
        @printf(io, "ORIGIN %f %f %f\n", o...)
        @printf(io, "SPACING %f %f %f\n", d...)
        @printf(io, "POINT_DATA %d\n", prod(n))
        println(io, "SCALARS Value float")
        println(io, "LOOKUP_TABLE default")
        for iz = 1:n[3], iy = 1:n[2], ix = 1:n[1]
            write(io, hton(Float32(H[iz, iy, ix])))
        end
    end
end

# * * * * * * * * * * * *
# *    Green Library    *
# * * * * * * * * * * * *

"""
    glib_readhead(io::IO) -> (n=n, x=x, y=y, z=z, t=t, risetime=rt)
"""
function glib_readhead(io::IO)
    rt = read(io, Float32)
    n = zeros(Int32, 4)
    read!(io, n)
    x = zeros(Float32, n[1])
    y = zeros(Float32, n[2])
    z = zeros(Float32, n[3])
    t = zeros(Float32, n[4])
    read!(io, x)
    read!(io, y)
    read!(io, z)
    read!(io, t)
    return (n=n, x=x, y=y, z=z, t=t, risetime=rt)
end

"""
    glib_readall(io::IO) -> (x, y, z, t, H, rt)
"""
function glib_readall(io::IO)
    (_, x, y, z, t, rt) = glib_readhead(io)
    H = zeros(Float32, length(t), 6, 3, length(z), length(y), length(x))
    read!(io, H)
    return (x, y, z, t, H, rt)
end

"""
    glib_readall(s::AbstractString) -> (x, y, z, t, H, rt)
"""
function glib_readall(s::AbstractString)
    return open(glib_readall, s, "r")
end

"""
    glib_readlocation(filename, x, y, z) -> (rt, t, w)
"""
function glib_readlocation(filename::AbstractString, x::Real, y::Real, z::Real)
    (n, xs, ys, zs, t, rt) = open(glib_readhead, filename, "r")
    if (x > maximum(xs)) || (x < minimum(xs)) || (y > maximum(ys)) || (y < minimum(ys)) || (z > maximum(zs)) ||
       (z < minimum(zs))
        error("Locaion out of range, require x($(minimum(xs)),$(maximum(xs))), y($(minimum(ys)),$(maximum(ys))), \
            z($(minimum(zs)),$(maximum(zs))), current is x:$x, y:$y, z:$z")
    end
    ix = max(2, findfirst(>(x), xs))
    iy = max(2, findfirst(>(y), ys))
    iz = max(2, findfirst(>(z), zs))
    h = (x - xs[ix-1]) / (xs[ix] - xs[ix-1])
    k = (y - ys[iy-1]) / (ys[iy] - ys[iy-1])
    l = (z - zs[iz-1]) / (zs[iz] - zs[iz-1])
    w = zeros(Float32, Int(n[4]), 6, 3)
    io = open(filename, "r")
    H = Mmap.mmap(io, Array{Float32,6}, (Int(n[4]), 6, 3, Int(n[3]), Int(n[2]), Int(n[1])),
        Int((sum(n) + 5) * 4))
    # for id = 1:3, ic = 1:6, it = 1:Int(n[4])
    #     w[it, ic, id] = H[it, ic, id, iz, iy, ix]
    # end
    for id = 1:3, ic = 1:6, it = 1:Int(n[4])
        w[it, ic, id] = H[it, ic, id, iz, iy, ix] * h * k * l +
                        H[it, ic, id, iz, iy, ix-1] * (1.0 - h) * k * l +
                        H[it, ic, id, iz, iy-1, ix] * h * (1.0 - k) * l +
                        H[it, ic, id, iz, iy-1, ix-1] * (1.0 - h) * (1.0 - k) * l +
                        H[it, ic, id, iz-1, iy, ix] * h * k * (1.0 - l) +
                        H[it, ic, id, iz-1, iy, ix-1] * (1.0 - h) * k * (1.0 - l) +
                        H[it, ic, id, iz-1, iy-1, ix] * h * (1.0 - k) * (1.0 - l) +
                        H[it, ic, id, iz-1, iy-1, ix-1] * (1.0 - h) * (1.0 - k) * (1.0 - l)
    end
    close(io)
    return (rt, t, w)
end

function glib_slice!(G::AbstractArray, xr::Tuple{Int,Int}, yr::Tuple{Int,Int}, zr::Tuple{Int,Int}, H::AbstractArray)
    nx = xr[2] - xr[1] + 1
    ny = yr[2] - yr[1] + 1
    nz = zr[2] - zr[1] + 1
    size(G) == (size(H, 1), size(H, 2), size(H, 3), nz, ny, nx)
    for ix = 1:nx, iy = 1:ny, iz = 1:nz
        G[:, :, :, iz, iy, ix] .= H[:, :, :, iz-1+zr[1], iy-1+yr[1], ix-1+xr[1]]
    end
    return nothing
end

function glib_writeglib(filename::AbstractString, H::AbstractArray{Float32},
    x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
    z::AbstractVector{<:Real}, t::AbstractVector{<:Real},
    rt::Real)
    if ndims(H) != 6
        error("dimenson of H not correct")
    end
    if size(H, 1) != length(t)
        error("size of t mismatch")
    end
    if size(H, 2) != 6
        error("size of component mismatch")
    end
    if size(H, 3) != 3
        error("size of channel mismatch")
    end
    if size(H, 4) != length(z)
        error("size of z mismatch")
    end
    if size(H, 5) != length(y)
        error("size of y mismatch")
    end
    if size(H, 6) != length(x)
        error("size of x mismatch")
    end
    open(filename, "w") do io
        write(io, Float32(rt))
        write(io, Int32(length(x)))
        write(io, Int32(length(y)))
        write(io, Int32(length(z)))
        write(io, Int32(length(t)))
        write(io, Float32.(x))
        write(io, Float32.(y))
        write(io, Float32.(z))
        write(io, Float32.(t))
        write(io, Float32.(H))
    end
end

function _cglib_page_size(b::AbstractVector{<:Integer})
    p = zeros(Int, length(b))
    for i = eachindex(b)
        if i == 1
            p[i] = 1
        else
            p[i] = p[i-1] * b[i-1]
        end
    end
    return p
end

function _cglib_lin2cart(l::Int, b::AbstractVector{<:Integer})
    c = zeros(Int, length(b))
    p = _cglib_page_size(b)
    i = length(b)
    res = l - 1
    while i > 0
        (d, r) = divrem(res, p[i])
        c[i] = Int(d)
        res = r
        i -= 1
    end
    return c .+ 1
end

@inline function cglib_predict(x1::Real, x2::Real, x3::Real, c::Integer)
    if c == 0
        return Float64(x3)
    elseif c == 1
        return 2 * Float64(x3) - Float64(x2)
    elseif c == 2
        return Float64(x1) - 3 * Float64(x2) + 3 * Float64(x3)
    else
        return Float64(x2)
    end
end

@inline function cglib_decodenum(r::Unsigned, expmask::Unsigned,
    nsig::Integer, sigmask::Unsigned, expshift::Integer)
    a = r & sigmask
    b = (r >> nsig) & expmask
    e = Int(b) - expshift
    if a < 2^(nsig - 1)
        s = Float64(a) * 2.0^(e - nsig + 1)
    else
        s = (Float64(a) - 2.0^(nsig)) * 2.0^(e - nsig + 1)
    end
    return s
end

"""
cglib_readhead(io::IO) ->
(nx, ny, nz, x0, y0, z0, dx, dy, dz, nt, dt, stf,
grid, 
bit_exp, expmask, expshift, bit_sig, sigmask, btype,
nleaf, leafv, nnode, leftnodes, rightnodes)
"""
function cglib_readhead(io::IO)
    flag = read(io, UInt8)
    partstart = zeros(UInt64, 4); read!(io, partstart)
    seek(io, 1 + 8 * 4 + 3 * 4)
    nt = read(io, Int32)
    dt = read(io, Float32)
    stf = zeros(Float32, nt); read!(io, stf)
    seek(io, partstart[1])
    ns = zeros(Int32, 3); read!(io, ns)
    xs = zeros(Float32, 6); read!(io, xs)
    (nx, ny, nz) = ns
    gridsize = Int.((nz, ny, nx))
    tracepospos = read(io, UInt64)
    seek(io, tracepospos)
    grid = zeros(UInt64, nz, ny, nx); read!(io, grid);

    seek(io, partstart[2])
    bit_exp = read(io, Int8)
    bit_sig = read(io, Int8)
    expshift = read(io, Int32)
    shiftval = read(io, Float64)
    nbit = bit_exp + bit_sig + 2

    typelist = (UInt8, UInt16, UInt32, UInt64)
    nbyte = ceil(Int, nextpow(2, nbit) / 8)
    btype = typelist[round(Int, log2(nbyte))+1]
    expmask = btype(2^bit_exp - 1)
    sigmask = btype(2^bit_sig - 1)

    nleaf = read(io, Int32)
    leafv = zeros(btype, nleaf); read!(io, leafv);
    nnode = read(io, Int32)
    leftnodes = zeros(Int32, nnode); read!(io, leftnodes);
    rightnodes = zeros(Int32, nnode); read!(io, rightnodes);

    return (nx, ny, nz, xs..., nt, dt, stf, grid, 
        bit_exp, expmask, expshift, bit_sig, sigmask, btype,
        nleaf, leafv, nnode, leftnodes, rightnodes)
end

const cglib_bitorflag = (0b10000000,
             0b01000000,
             0b00100000,
             0b00010000,
             0b00001000,
             0b00000100,
             0b00000010,
             0b00000001);

function cglib_readtrace(io::IO, ix::Integer, iy::Integer, iz::Integer,
    nt::Integer, grid::Array{UInt64,3},
    nnode::Integer, leftnodes, rightnodes, leafv,
    bit_exp::Integer, expmask::Unsigned, expshift::Integer, bit_sig::Integer, sigmask::Unsigned
    )
    seek(io, grid[iz,iy,ix])
    tp = read(io, Float32)
    ts = read(io, Float32)
    nzero = read(io, Int32)
    amp = read(io, Float32)
    cbyte = read(io, Int32)
    cbits = read(io, Int8)
    encoded = zeros(UInt8, cbyte); read!(io, encoded);
    
    decoded = zeros(Float32, nt, 6, 3);
    tracedatasize = [nt-nzero, 6, 3]
    idata = 1
    ibyte = 0
    ibit = 8
    inode = nnode
    while true
        ibit += 1
        if ibit == 9
            ibit = 1
            ibyte += 1
        end
        inode = iszero(encoded[ibyte] & cglib_bitorflag[ibit]) ? leftnodes[inode] : rightnodes[inode]
        if iszero(leftnodes[inode])
            # println(idata)
            (it, im, ic) = _cglib_lin2cart(idata, tracedatasize)
            it += nzero
            pretype = (leafv[inode] >> (bit_exp+bit_sig)) & 0b11;
            if it == 1
                hp = 0.0
            elseif it == 2
                hp = cglib_predict(0.0, 0.0, decoded[it-1,im,ic], pretype)
            elseif it == 3
                hp = cglib_predict(0.0, decoded[it-2,im,ic], decoded[it-1,im,ic], pretype)
            else
                hp = cglib_predict(decoded[it-3,im,ic], decoded[it-2,im,ic], decoded[it-1,im,ic], pretype)
            end
            # decoded[it, im, ic] =  hp + dr * leafv[inode] - shiftval
            decoded[it,im,ic] = hp + cglib_decodenum(leafv[inode], expmask, bit_sig, sigmask, expshift)*Float64(amp);
            idata += 1
            inode = nnode
        end
        if ibyte == cbyte && ibit == cbits
            break
        end
    end
    return (tp, ts, decoded);
end

"""
cglib_readlocation(filename, x, y, z) -> (stf, dt, tp, ts, g)
"""
function cglib_readlocation(filename::AbstractString, x::Real, y::Real, z::Real)
    io = open(filename);
    (nx, ny, nz, x0, y0, z0, dx, dy, dz, nt, dt, stf, grid, 
        bit_exp, expmask, expshift, bit_sig, sigmask, btype,
        nleaf, leafv, nnode, leftnodes, rightnodes) = cglib_readhead(io);
    if (x > (x0+(nx-1)*dx)) || (x < x0) ||
        (y > (y0+(ny-1)*dy)) || (y < y0) ||
        (z > (z0+(nz-1)*dz)) || (z < z0)
         error("Locaion out of range, require x($(x0),$(x0+(nx-1)*dx)), y($(y0),$(y0+(ny-1)*dy)), \
             z($(z0),$(z0+(nz-1)*dz)), current is x:$x, y:$y, z:$z")
     end
    xp = floor(Int, (x - x0) / dx) + 1;
    yp = floor(Int, (y - y0) / dy) + 1;
    zp = floor(Int, (z - z0) / dz) + 1;
    if xp == nx
        xp = nx-1;
    end
    if yp == ny
        yp = ny-1;
    end
    if zp == nz
        zp = nz - 1;
    end
    h = (x - x0) / dx - xp + 1.0
    k = (y - y0) / dy - yp + 1.0
    l = (z - z0) / dz - zp + 1.0
    parp = (nt, grid, nnode, leftnodes, rightnodes, leafv, bit_exp, expmask, expshift, bit_sig, sigmask)
    w = zeros(Float32, nt, 6, 3)
    gt1 = cglib_readtrace(io, xp,   yp,   zp,   parp...)
    gt2 = cglib_readtrace(io, xp+1, yp,   zp,   parp...)
    gt3 = cglib_readtrace(io, xp,   yp+1, zp,   parp...)
    gt4 = cglib_readtrace(io, xp+1, yp+1, zp,   parp...)
    gt5 = cglib_readtrace(io, xp,   yp,   zp+1, parp...)
    gt6 = cglib_readtrace(io, xp+1, yp,   zp+1, parp...)
    gt7 = cglib_readtrace(io, xp,   yp+1, zp+1, parp...)
    gt8 = cglib_readtrace(io, xp+1, yp+1, zp+1, parp...)
    w .+= gt1[3] .* ((1.0-h) * (1.0-k) * (1.0-l))
    w .+= gt2[3] .* (h       * (1.0-k) * (1.0-l))
    w .+= gt3[3] .* ((1.0-h) * k       * (1.0-l))
    w .+= gt4[3] .* (h       * k       * (1.0-l))
    w .+= gt5[3] .* ((1.0-h) * (1.0-k) * l)
    w .+= gt6[3] .* (h       * (1.0-k) * l)
    w .+= gt7[3] .* ((1.0-h) * k       * l)
    w .+= gt8[3] .* (h       * k       * l)
    tp = gt1[1] * ((1.0-h) * (1.0-k) * (1.0-l)) +
         gt2[1] * (h       * (1.0-k) * (1.0-l)) +
         gt3[1] * ((1.0-h) * k       * (1.0-l)) +
         gt4[1] * (h       * k       * (1.0-l)) +
         gt5[1] * ((1.0-h) * (1.0-k) * l) +
         gt6[1] * (h       * (1.0-k) * l) +
         gt7[1] * ((1.0-h) * k       * l) +
         gt8[1] * (h       * k       * l)
    ts = gt1[2] * ((1.0-h) * (1.0-k) * (1.0-l)) +
        gt2[2] * (h       * (1.0-k) * (1.0-l)) +
        gt3[2] * ((1.0-h) * k       * (1.0-l)) +
        gt4[2] * (h       * k       * (1.0-l)) +
        gt5[2] * ((1.0-h) * (1.0-k) * l) +
        gt6[2] * (h       * (1.0-k) * l) +
        gt7[2] * ((1.0-h) * k       * l) +
        gt8[2] * (h       * k       * l)
    close(io)
    return (stf, dt, tp, ts, w);
end

# * * * * * * * * * * * *
# *        Model        *
# * * * * * * * * * * * *

"""
"""
function model_writebin(fn::AbstractString, vp::AbstractArray, vs::AbstractArray, rho::AbstractArray,
    topz::Real, dx::Real, dy::Real, dz::Real, airmask::AbstractArray{Bool}=Bool[])
    (v0, dv, nsp, materialhash, materialtable) = discretematerial((vp, vs, rho), (100, 100, 20), .!airmask)
    materialgrid = field2material((vp, vs, rho), v0, dv, materialhash, airmask)
    open(fn, "w") do io
        # * field
        write(io, Int32(length(v0)))
        write(io, Float32.(v0))
        write(io, Float32.(dv))
        write(io, Int32.(nsp))
        write(io, Int32(size(materialtable, 2)))
        write(io, Float32.(materialtable))
        # * mesh
        write(io, Int32(size(vp, 3)))
        write(io, Int32(size(vp, 2)))
        write(io, Int32(size(vp, 1)))
        write(io, Float32(dx))
        write(io, Float32(dy))
        write(io, Float32(dz))
        write(io, Float32(topz))
        # * model
        write(io, Int32.(materialgrid))
    end
    return nothing
end
