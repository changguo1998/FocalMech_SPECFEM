using Printf, LinearAlgebra

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

"""
    ijk2n(ijk::AbstractVector, base::Tuple)

Cartesian Index -> Linear Index
"""
function ijk2n(ijk::AbstractVector, base::Tuple)
    n = 1
    b = 1
    for i in eachindex(base)
        n += (ijk[i] - 1) * b
        b *= base[i]
    end
    return Int(n)
end

"""
    n2ijk(n::Int, base::Tuple)

Linear Index -> Cartesian Index
"""
function n2ijk(n::Int, base::Tuple)
    t = zeros(Int, length(base))
    bcum = 1
    scum = 1
    for i in eachindex(base)
        t[i] = mod((n - scum) / bcum, base[i]) + 1 |> Int
        scum += (t[i] - 1) * bcum
        bcum *= base[i]
    end
    return Tuple(t)
end

"""
    frameinterp!(v::AbstractArray, mask::AbstractArray{Bool}, ndim::Int)

fill unmasked value using linear interpolation
"""
function frameinterp!(v::AbstractArray, mask::AbstractArray{Bool}, ndim::Int)
    vsize = size(v)
    othermesh = vsize[[1:ndim-1; ndim+1:end]]
    dimlen = vsize[ndim]
    Threads.@threads for idx in CartesianIndices(othermesh)
        maskidx = map(i -> CartesianIndex(idx.I[1:ndim-1]..., i, idx.I[ndim:end]...), 1:dimlen)
        p = findnext(mask[maskidx], 1)
        q = findnext(mask[maskidx], p + 1)
        while !isnothing(q)
            l = q - p
            for s = 1:l-1
                v[idx.I[1:ndim-1]..., p+s, idx.I[ndim:end]...] = v[idx.I[1:ndim-1]..., p, idx.I[ndim:end]...] +
                                                                 (v[idx.I[1:ndim-1]..., q, idx.I[ndim:end]...] -
                                                                  v[idx.I[1:ndim-1]..., p, idx.I[ndim:end]...]) * s / l
            end
            p = q
            q = findnext(mask[maskidx], p + 1)
        end
    end
    return nothing
end

function smooth!(y::AbstractMatrix, x::AbstractMatrix, r::Int=2)
    (m, n) = size(x)
    # for i = 1:m, j = 1:n
    Threads.@threads for idx = CartesianIndices(x)
        (i, j) = idx.I
        c = 0.0
        y[i, j] = 0.0
        for p = -r:r, q = -r:r
            if (p + i >= 1) && (p + i <= m) && (q + j >= 1) && (q + j <= n)
                c += 1.0
                y[i, j] += x[i+p, j+q]
            end
        end
        y[i, j] /= c
    end
    return nothing
end

function smoothmeshcoor!(mesh::AbstractArray, dim::Int=1)
    msize = size(mesh)
    othermesh = msize[[1:dim-1; dim+1:end]]
    buf = zeros(othermesh...)
    for i = 1:size(mesh, dim)
        n = selectdim(mesh, dim, i)
        buf .= n
        smooth!(n, buf)
    end
    return nothing
end

"""
    geninterval(nint::Int, v::AbstractArray)

generate parameter to split field into several intervals
nint: number of interval
v:    values
"""
function geninterval(nint::Int, v::AbstractArray)
    vdiff = maximum(v) - minimum(v)
    if vdiff > 0.0
        if nint > 1
            dv = (maximum(v) - minimum(v)) / (nint - 1)
        else
            dv = (maximum(v) - minimum(v)) * 1.2
        end
        v0 = minimum(v) - dv / 2
    else
        dv = 1.0
        v0 = minimum(v) - dv / 2
    end
    return (v0, dv, floor(Int, (maximum(v) - v0) / dv) + 1)
end

"""
    discretematerial(fields::Tuple, nints::Tuple, mask::AbstractArray{Bool}=Bool[])

fields  different field data like vp, vs
nints   number of intervals for each field
mask    datas that will be used in calculation
"""
function discretematerial(fields::Tuple, nints::Tuple, mask::AbstractArray{Bool}=Bool[])
    nfield = length(nints)
    v0 = zeros(nfield)
    dv = zeros(nfield)
    nsplit = zeros(Int, nfield)
    if isempty(mask)
        for i in eachindex(nints)
            (v0[i], dv[i], nsplit[i]) = geninterval(nints[i], fields[i])
        end
    else
        for i in eachindex(nints)
            (v0[i], dv[i], nsplit[i]) = geninterval(nints[i], fields[i][mask])
        end
    end
    materialhash = zeros(Int, nsplit...)
    for i in eachindex(materialhash)
        materialhash[i] = i
    end
    mattable = zeros(nfield, length(materialhash))
    for idx in CartesianIndices(Tuple(nsplit))
        tableid = materialhash[idx]
        for i = 1:nfield
            mattable[i, tableid] = v0[i] + (idx.I[i] - 0.5) * dv[i]
        end
    end
    return (v0, dv, nsplit, materialhash, mattable)
end

"""
    field2material(fields, v0, dv, materialhash, mask)

transformat field data to material number. when mask is `true`, the point will be filled with 0
"""
function field2material(fields::Tuple, v0::AbstractVector, dv::AbstractVector, materialhash::AbstractArray,
                        mask::AbstractArray{Bool}=Bool[])
    nfield = length(v0)
    materialfield = zeros(Int, size(fields[1]))
    # println("v0: ", v0)
    # println("dv: ", dv)
    usemask = !isempty(mask)
    Threads.@threads for idx in eachindex(materialfield)
        if usemask && mask[idx]
            materialfield[idx] = 0
        else
            t = zeros(Int, nfield)
            for i = 1:nfield
                # if isnan(fields[i][idx])
                #     println(i, " ", idx, " ", fields[i][idx])
                #     error("2")
                # end
                t[i] = floor(Int, (fields[i][idx] - v0[i]) / dv[i]) + 1
                # if t[i] < 1
                #     println("i:     ", i)
                #     println("idx:   ", idx)
                #     println("field: ", fields[i][idx])
                #     println("v0:    ", v0[i])
                #     println("dv:    ", dv[i])
                #     error("index <= 0")
                # end
            end
            materialfield[idx] = materialhash[t[1], t[2], t[3]]
        end
    end
    return materialfield
end

"""
    maskgridinterface(mat::AbstractArray)

set the interface point to `true`
"""
function maskgridinterface(mat::AbstractArray)
    dimsize = size(mat)
    mask = falses(dimsize)
    Threads.@threads for idx in CartesianIndices(dimsize[2:end])
        mask[1, idx] = true
        mask[end, idx] = true
        for iz = 2:dimsize[1]-1
            mask[iz, idx] = mat[iz, idx] != mat[iz-1, idx]
        end
    end
    return mask
end

@inline transindex(i, b1, b2) = round(Int, (i - 1) / (b1 - 1) * (b2 - 1)) + 1

"""
    shiftmesh_alldirc!(mask, maskgrid)
"""
function shiftmesh_alldirc!(mask::AbstractArray{Bool}, maskgrid::AbstractArray{Bool})
    @inline dist(i, bi, j, bj) = sum(@.(((i - 1) / (bi - 1) - (j - 1) / (bj - 1))^2))
    ndim = ndims(mask)
    ssize = size(maskgrid)
    tsize = size(mask)
    binding = zeros(Int, ndim, size(mask)...)
    rmin = fill(-1.0, tsize)
    r0 = dist(fill(2.0, ndim), ssize, ones(ndim), ssize)
    cmax = zeros(Int, ndim)
    cmin = zeros(Int, ndim)
    for idx1 in CartesianIndices(maskgrid)
        if maskgrid[idx1]
            for i = 1:ndim
                t = transindex(idx1.I[i], ssize[i], tsize[i])
                cmin[i] = max(1, t - 1)
                cmax[i] = min(tsize[i], t + 1)
            end
            for idx2 in CartesianIndices(Tuple(map((p, q) -> p:q, cmin, cmax)))
                r = dist(idx1.I, ssize, idx2.I, tsize)
                if (r < r0) && ((r < rmin[idx2]) || (rmin[idx2] < 0.0))
                    for i = 1:ndim
                        binding[i, idx2] = idx1.I[i]
                    end
                    rmin[idx2] = r
                end
            end
        end
    end
    mesh = zeros(ndim, tsize...)
    for idx in CartesianIndices(mask)
        if rmin[idx] >= 0.0
            mask[idx] = true
            for i = 1:ndim
                mesh[i, idx] = (binding[i, idx] - 1) / (ssize[i] - 1)
            end
        else
            mask[idx] = false
        end
    end
    return mesh
end

"""
shiftmesh_zdirc!(mask, maskgrid)
"""
function shiftmesh_zdirc!(mask::AbstractArray{Bool}, maskgrid::AbstractArray{Bool})
    ndim = ndims(mask)
    ssize = size(maskgrid)
    tsize = size(mask)
    mesh = zeros(ndim, tsize...)
    mask .= false
    shiftlim = max(1, round(Int, (size(maskgrid, 1) - 1) / (size(mask, 1) - 1) / 3))
    Threads.@threads for idx1 in CartesianIndices(tsize[2:end])
        tidx = zeros(Int, ndim - 1)
        for i = 1:ndim-1
            tidx[i] = transindex(idx1.I[i], tsize[i+1], ssize[i+1])
        end
        idx2 = CartesianIndex(tidx...)
        for jz = 1:ssize[1]
            if maskgrid[jz, idx2]
                iz = transindex(jz, ssize[1], tsize[1])
                if abs(transindex(iz, tsize[1], ssize[1]) - jz) <= shiftlim
                    mesh[1, iz, idx1] = (jz - 1) / (ssize[1] - 1)
                else
                    mesh[1, iz, idx1] = (iz - 1) / (tsize[1] - 1)
                end
                for i = 2:ndim
                    mesh[i, iz, idx1] = (idx1.I[i-1] - 1) / (tsize[i] - 1)
                end
                mask[iz, idx1] = true
            end
        end
    end
    return mesh
end

"""
shiftmesh_onlytopo!(mask, maskgrid)
"""
function shiftmesh_onlytopo!(mask::AbstractArray{Bool}, maskgrid::AbstractArray{Bool})
    ndim = ndims(mask)
    ssize = size(maskgrid)
    tsize = size(mask)
    mesh = zeros(ndim, tsize...)
    mask .= false
    # Threads.@threads for idx in CartesianIndices(mask)
    #     if any(==(1), idx.I) || any(tsize .== idx.I)
    #         for d = 1:ndim
    #             mesh[d, idx] = (idx.I[d] - 1.0) / (tsize[d] - 1.0)
    #         end
    #         mask[idx] = true
    #     end
    # end
    Threads.@threads for idx in CartesianIndices(tsize[2:end])
        for d = 2:ndim
            mesh[d, 1, idx] = (idx.I[d-1] - 1.0) / (tsize[d] - 1.0)
            mesh[d, end, idx] = (idx.I[d-1] - 1.0) / (tsize[d] - 1.0)
        end
        mesh[1, 1, idx] = 0.0
        mesh[1, end, idx] = 1.0
        mask[1, idx] = true
        mask[end, idx] = true
    end
    Threads.@threads for idx1 in CartesianIndices(tsize[2:end])
        tidx = zeros(Int, ndim - 1)
        for i = 1:ndim-1
            tidx[i] = transindex(idx1.I[i], tsize[i+1], ssize[i+1])
        end
        idx2 = CartesianIndex(tidx...)
        count = 0
        jz = 0
        for jjz = 1:ssize[1]
            if maskgrid[jjz, idx2]
                count += 1
            end
            if count == 2
                jz = jjz
                break
            end
        end
        mesh[1, 2, idx1] = (jz - 1) / (ssize[1] - 1)
        for d = 2:ndim
            mesh[d, 2, idx1] = (idx1.I[d-1] - 1) / (tsize[d] - 1)
        end
        mask[2, idx1] = true
    end
    return mesh
end

"""
shiftmesh_split!(mask, maskgrid)
"""
function shiftmesh_split!(mask::AbstractArray{Bool}, maskgrid::AbstractArray{Bool})
    ndim = ndims(mask)
    ssize = size(maskgrid)
    tsize = size(mask)
    mesh = zeros(ndim, tsize...)
    mask .= false
    nlayer = size(mask, 1) - 1
    # shiftlim = max(1, round(Int, (size(maskgrid, 1)-1)/(size(mask, 1)-1)/2))
    Threads.@threads for idx1 in CartesianIndices(tsize[2:end])
        # target index
        tidx = zeros(Int, ndim - 1)
        for i = 1:ndim-1
            tidx[i] = transindex(idx1.I[i], tsize[i+1], ssize[i+1])
        end
        idx2 = CartesianIndex(tidx...)
        # depth of interface in source model
        interfacez = zeros(ssize[1])
        nintf = 0
        for jz = 1:ssize[1]
            if maskgrid[jz, idx2]
                nintf += 1
                interfacez[nintf] = (jz - 1) / (ssize[1] - 1)
            end
        end
        # interface depth -> layer thickness
        nthick = nintf - 1
        thickness = zeros(nthick)
        for i = 1:nthick
            thickness[i] = interfacez[i+1] - interfacez[i]
        end
        # merge redundant layers
        while nthick > nlayer
            # try each way of merge and find the thinnest one
            p = 0
            v = 0.0
            for q = 2:nthick-1
                if (thickness[q] + thickness[q+1] < v) || iszero(v)
                    p = q
                    v = thickness[q] + thickness[q+1]
                end
            end
            # merge thickness
            thickness[p] += thickness[p+1]
            for i = p+1:nthick-1
                thickness[i] = thickness[i+1]
            end
            thickness[nthick] = 0.0
            # merge depth
            for i = p+1:nintf-1
                interfacez[i] = interfacez[i+1]
            end
            interfacez[nintf] = 0.0
            nintf -= 1
            nthick -= 1
        end

        nelementz = floor.(Int, thickness .* nlayer)
        nelementz[1] = 1
        for i = 1:nthick
            if nelementz[i] < 1
                nelementz[i] = 1
            end
        end
        # add more layer when number of layer is less than maskgrid
        while sum(nelementz) < nlayer
            avethick = thickness ./ nelementz
            adjustable = trues(length(avethick))
            adjustable[1] = false
            p = 0
            v = 0.0
            for q = 2:nthick
                if adjustable[q] && ((avethick[q] > v) || iszero(v))
                    p = q
                    v = avethick[q]
                end
            end
            nelementz[p] += 1
        end
        # remove layer when number of layer is larger than maskgrid
        while sum(nelementz) > nlayer
            avethick = thickness ./ nelementz
            adjustable = nelementz .> 1
            adjustable[1] = false
            p = 0
            v = 0.0
            for q = 2:nthick
                if adjustable[q] && ((avethick[q] < v) || iszero(v))
                    p = q
                    v = avethick[q]
                end
            end
            nelementz[p] -= 1
        end
        p = 1
        for i = 1:nthick
            mesh[1, p, idx1] = interfacez[i]
            for j = 2:ndim
                mesh[j, p, idx1] = (idx1.I[j-1] - 1) / (tsize[j] - 1)
            end
            mask[p, idx1] = true
            p += nelementz[i]
        end
        # @info "sum: $(sum(nelementz)), final coor: $(p), final point: $(nintf)"
        mesh[1, p, idx1] = interfacez[nintf]
        for i = 2:ndim
            mesh[i, p, idx1] = (idx1.I[i-1] - 1) / (tsize[i] - 1)
        end
        mask[p, idx1] = true
    end
    return mesh
end

function fillmaterial_select(semmesh::AbstractArray, material_grid::AbstractArray)
    t = size(semmesh)
    ssize = size(material_grid)
    tsize = t[2:end]
    material = zeros(Int, tsize)
    ndim = t[1]
    idx2 = zeros(Int, ndim)
    for idx in CartesianIndices(tsize)
        for i = 1:ndim
            idx2[i] = max(round(Int, semmesh[i, idx] * (ssize[i] - 1)) + 1, 1)
        end
        material[idx] = material_grid[idx2...]
    end
    return material
end

function dualvector!(y::AbstractMatrix, x::AbstractMatrix)
    y[:, 1] = cross(x[:, 2], x[:, 3])
    y[:, 2] = cross(x[:, 3], x[:, 1])
    y[:, 3] = cross(x[:, 1], x[:, 2])
    for i = 1:3
        if dot(x[:, i], y[:, i]) < 0.0
            y[:, i] .*= -1.0
        end
        normalize!(@view(y[:, i]))
    end
    return nothing
end

@inline δ(i, j) = (i == j) ? 1 : 0

function interpfield2cell(semmesh::AbstractArray, fields::Tuple)
    t = size(semmesh)
    ndim = t[1]
    tsize = t[2:end]
    csize = tsize .- 1
    nfield = length(fields)
    ssize = size(fields[1])
    cellmat = Vector{Array{Float64,ndim}}(undef, nfield)
    for iv = 1:nfield
        cellmat[iv] = zeros(csize)
    end
    Threads.@threads for idx1 in CartesianIndices(csize)
        tvec = zeros(ndim, ndim)
        dual1 = zeros(ndim, ndim)
        dual2 = zeros(ndim, ndim)
        tidx1 = zeros(Int, ndim)
        tidx2 = zeros(Int, ndim)
        tid = zeros(Int, ndim)
        cmin = zeros(Int, ndim)
        cmax = zeros(Int, ndim)
        vcums = zeros(nfield)
        # cmin .= 0
        # cmax .= 0
        for i = 1:ndim
            tidx1[i] = max(round(Int, semmesh[i, idx1] * (ssize[i] - 1)) + 1, 1)
            if (tidx1[i] < cmin[i]) || (cmin[i] < 1)
                cmin[i] = tidx1[i]
            end
            if (tidx1[i] > cmax[i]) || (cmax[i] < 1)
                cmax[i] = tidx1[i]
            end
        end
        for i = 1:ndim
            for j = 1:ndim
                tid[j] = idx1.I[j] + δ(i, j)
            end
            for j = 1:ndim
                ti = max(1, round(Int, semmesh[j, tid...] * (ssize[j] - 1)) + 1)
                tvec[j, i] = ti - tidx1[j]
                if (ti < cmin[j]) || (cmin[j] < 1)
                    cmin[j] = ti
                end
                if (ti > cmax[j]) || (cmax[j] < 1)
                    cmax[j] = ti
                end
            end
        end
        dualvector!(dual1, tvec)
        for i = 1:ndim
            tid[i] = idx1.I[i] + 1
        end
        for i = 1:ndim
            tidx2[i] = max(1, round(Int, semmesh[i, tid...] * (ssize[i] - 1)) + 1)
            if (tidx2[i] < cmin[i]) || (cmin[i] < 1)
                cmin[i] = tidx2[i]
            end
            if (tidx2[i] > cmax[i]) || (cmax[i] < 1)
                cmax[i] = tidx2[i]
            end
        end
        for i = 1:ndim
            for j = 1:ndim
                tid[j] = idx1.I[j] - δ(i, j) + 1
            end
            for j = 1:ndim
                ti = max(1, round(Int, semmesh[j, tid...] * (ssize[j] - 1)) + 1)
                tvec[j, i] = ti - tidx2[j]
                if (ti < cmin[j]) || (cmin[j] < 1)
                    cmin[j] = ti
                end
                if (ti > cmax[j]) || (cmax[j] < 1)
                    cmax[j] = ti
                end
            end
        end
        dualvector!(dual2, tvec)
        vcums .= 0.0
        ncums = 0.0
        for idx2 in CartesianIndices(Tuple(map((p, q) -> p:q, cmin, cmax)))
            flag = true
            for i = 1:3
                flag &= dot(dual1[:, i], idx2.I .- tidx1) >= 0.0
                flag &= dot(dual2[:, i], idx2.I .- tidx1) >= 0.0
            end
            if flag
                for iv = 1:nfield
                    vcums[iv] += fields[iv][idx2]
                end
                ncums += 1.0
            end
        end
        for iv = 1:nfield
            cellmat[iv][idx1] = vcums[iv] / max(1.0, ncums)
        end
    end
    return Tuple(cellmat)
end

function locatelin(x::Real, xs::Vector{<:Real})
    if (x < xs[1]) || (x > xs[end])
        println("x: ", x, " xs: ", xs)
        error("out of range")
    end
    # i = findlast(<(x), xs)
    i = findfirst(>(x), xs)
    if isnothing(i)
        i = length(xs) - 1
    else
        i -= 1
    end
    return (i, (x-xs[i])/(xs[i+1]-xs[i]))
end

function locatepoint(coor, mesh, x, y)
    # println(coor)
    msize = size(mesh)
    id = zeros(Int, 3)
    (p, h) = locatelin(coor[3], x)
    (q, k) = locatelin(coor[2], y)
    id[3] = p
    id[2] = q
    z = zeros(msize[2])
    for i = 1:msize[2]
        # tz = 0.0
        # for (m, n) in ([0, 0], [0, 1], [1, 0], [1, 1])
        #     tz += mesh[1, i, id[2]+m, id[3]+n] * p[1, m+1] * p[2, n+1]
        # end
        z[i] = mesh[1, i, q, p] * (1.0 - k) * (1.0 - h) +
             mesh[1, i, q+1, p] * k * (1.0 - h) +
             mesh[1, i, q, p+1] * (1.0 - k) * h +
             mesh[1, i, q+1, p+1] * k * h
    end
    z[1] = 0.0
    z[end] = 1.0
    # iz -= 1
    # id[1] = iz
    (id[1], r) = locatelin(coor[1], z)
    # println(id)
    return CartesianIndex(id[1], id[2], id[3])
end

function _nearestpoint(coor, mesh, fdsize)
    msize = size(mesh)
    if coor.I[3] == msize[4]
        mx = coor.I[3]
    else
        mx = coor.I[3] + 0.5
    end
    if coor.I[2] == msize[3]
        my = coor.I[2]
    else
        my = coor.I[2] + 0.5
    end
    mz = 0.0
    cz = 0
    for sx = 0:1, sy = 0:1, sz = 0:1
        if (coor.I[1] == msize[2]) && (sz == 1)
            continue
        end
        if (coor.I[2] == msize[3]) && (sy == 1)
            continue
        end
        if (coor.I[3] == msize[4]) && (sx == 1)
            continue
        end
        mz += mesh[1, coor.I[1]+sz, coor.I[2]+sy, coor.I[3]+sx]
        cz += 1
    end
    mz /= cz
    return CartesianIndex(round(Int, mz*fdsize[1]), round(Int, (my-1)/(msize[3]-1)*(fdsize[2]-1))+1,
        round(Int, (mx-1)/(msize[4]-1)*(fdsize[3]-1))+1)
end

function interpfield2cell_looponfd(semmesh::AbstractArray, fields::Tuple, fieldair::AbstractArray{<:Bool, 3})
    t = size(semmesh)
    ndim = t[1]
    tsize = t[2:end]
    csize = tsize .- 1
    nfield = length(fields)
    ssize = size(fields[1])
    cellmat = Vector{Array{Float64,ndim}}(undef, nfield)
    cellcount = zeros(csize)
    locks = Array{Threads.SpinLock, 3}(undef, csize)
    for idx in CartesianIndices(csize)
        locks[idx] = Threads.SpinLock()
    end
    for iv = 1:nfield
        cellmat[iv] = zeros(csize)
    end
    x = semmesh[3,1,1,:]
    y = semmesh[2,1,:,1]
    # println(x)
    # println(y)
    Threads.@threads for idx in CartesianIndices(ssize)
        if !fieldair[idx]
            normcoor = zeros(ndim)
            for i = 1:ndim
                normcoor[i] = (idx.I[i] - 1) / (ssize[i] - 1)
            end
            idx2 = locatepoint(normcoor, semmesh, x, y)
            # println(normcoor, " -> ", idx2)
            lock(locks[idx2]) do
                for iv = 1:nfield
                    cellmat[iv][idx2] += fields[iv][idx]
                end
                cellcount[idx2] += 1.0
            end
        end
    end
    GC.gc()
    # if any(iszero, cellcount[2:end,:,:])
    #     for idx = CartesianIndices(csize)
    #         if (idx.I[1] > 1) && iszero(cellcount[idx])
    #             println(idx)
    #         end
    #     end
    #     error("some cell not filled")
    # end
    Threads.@threads for idx in CartesianIndices(csize)
        if idx.I[1] == 1
            continue
        end
        if iszero(cellcount[idx])
            idx2 = _nearestpoint(idx, semmesh, ssize)
            for iv = 1:nfield
                cellmat[iv][idx] = fields[iv][idx2]
            end
            cellcount[idx] = 1.0
        end
    end
    GC.gc()
    Threads.@threads for idx in eachindex(cellcount)
        if iszero(cellcount[idx])
            cellcount[idx] = 1.0
        end
    end
    GC.gc()
    for iv = 1:nfield
        Threads.@threads for idx in eachindex(cellcount)
            cellmat[iv][idx] /= cellcount[idx]
        end
    end
    return Tuple(cellmat)
end

function interpunmaskedcoor!(mesh::AbstractArray, mask::AbstractArray{Bool})
    ndim = ndims(mask)
    nsize = size(mask)
    tmask = falses(nsize...)

    for i = 1:ndim
        tmask .= mask
        Threads.@threads for idx in CartesianIndices(tmask)
            if !tmask[idx]
                if idx.I[i] == 1
                    tmask[idx] = true
                    mesh[i, idx] = 0.0
                end
                if idx.I[i] == nsize[i]
                    tmask[idx] = true
                    mesh[i, idx] = 1.0
                end
            end
        end
        frameinterp!(selectdim(mesh, 1, i), tmask, i)
    end
    return nothing
end

function _cell2region(cell::AbstractArray)
    ndim = ndims(cell)
    csize = size(cell)
    flag = trues(csize)
    indexs = zeros(Int, csize)
    for i in eachindex(cell)
        indexs[i] = i
    end
    table    = zeros(Int, ndim * 2 + 1, 0)
    trange   = zeros(Int, ndim * 2 + 1)
    pertdirc = 0:2*ndim-1
    cmin     = zeros(Int, ndim)
    cmax     = zeros(Int, ndim)
    while any(flag)
        seed = rand(indexs[flag])
        seedijk = n2ijk(seed, csize)
        for i = 1:ndim
            trange[i*2-1] = seedijk[i]
            trange[i*2]   = seedijk[i]
        end
        trange[ndim*2+1] = cell[seed]
        changeflag = true
        while changeflag
            changeflag = false
            for i in pertdirc
                pertdim = mod(i, ndim) + 1
                pertdrc = (i < ndim) ? 1 : 0
                rid = pertdim * 2 - 1 + pertdrc
                rsft = pertdrc * 2 - 1
                for j = 1:ndim
                    if j == pertdim
                        cmin[j] = cmax[j] = max(1, min(csize[pertdim], trange[rid] + rsft))
                    else
                        cmin[j] = trange[j*2-1]
                        cmax[j] = trange[j*2]
                    end
                end
                flagnewcell = cmin[pertdim] != trange[rid]
                if !flagnewcell
                    continue
                end
                flagallsame = true
                for idx in CartesianIndices(Tuple(map((p, q) -> p:q, cmin, cmax)))
                    flagallsame &= isequal(cell[idx], trange[end])
                    if !flagallsame
                        break
                    end
                end
                if flagallsame && flagnewcell
                    trange[rid] = max(1, min(csize[pertdim], trange[rid] + rsft))
                    changeflag = true
                    break
                end
            end # for i = pertdirc
        end # while changeflag
        table = hcat(table, trange)
        for i = 1:ndim
            cmin[i] = trange[i*2-1]
            cmax[i] = trange[i*2]
        end
        for idx in CartesianIndices(Tuple(map((p, q) -> p:q, cmin, cmax)))
            flag[idx] = false
        end
    end # while any(flag)
    return table
end

function isprime(n::Int)
    t = floor(Int, sqrt(n))
    while t > 1
        if mod(n, t) == 0
            break
        end
        t -= 1
    end
    return t == 1
end

function firstprime(n::Int)
    if isprime(n)
        return n
    end
    for i = 2:n
        if mod(n, i) == 0
            return i
        end
    end
end

function decomposeprime(n::Int)
    c = Int[]
    t = n
    while t != 1
        push!(c, firstprime(t))
        t = round(Int, t/c[end])
    end
    return c
end

function _balance_decompose(n::Int, d::Int)
    buf = ones(Int, d)
    cmps = decomposeprime(n)
    l = length(cmps)
    for i in eachindex(cmps)
        (_, j) = findmin(buf)
        buf[j] *= cmps[l-i+1]
    end
    return buf
end

function cell2region(cell::AbstractArray)
    nth = Threads.nthreads()
    if nth == 1
        return _cell2region(cell)
    end
    ndim = ndims(cell)
    nslice = _balance_decompose(nth, ndim)
    resultpart = Vector{Matrix{Int}}(undef, nth)
    dx = round(Int, size(cell, 3)/nslice[1])
    xl = zeros(Int, nslice[1] + 1)
    for i = 1:nslice[1]
        xl[i] = (i - 1) * dx + 1
    end
    xl[end] = size(cell, 3)

    dy = round(Int, size(cell, 2)/nslice[2])
    yl = zeros(Int, nslice[2] + 1)
    for i = 1:nslice[2]
        yl[i] = (i - 1) * dy + 1
    end
    yl[end] = size(cell, 2)

    dz = round(Int, size(cell, 1)/nslice[3])
    zl = zeros(Int, nslice[3] + 1)
    for i = 1:nslice[3]
        zl[i] = (i - 1) * dz + 1
    end
    zl[end] = size(cell, 1)

    limitpart = Tuple{Int,Int,Int,Int,Int,Int}[]
    for ix = 1:nslice[1], iy = 1:nslice[2], iz = 1:nslice[3]
        push!(limitpart, (xl[ix], xl[ix+1], yl[iy], yl[iy+1], zl[iz], zl[iz+1]))
    end

    Threads.@threads for i in eachindex(limitpart)
        lim = limitpart[i]
        tc = deepcopy(cell[lim[5]:lim[6], lim[3]:lim[4], lim[1]:lim[2]])
        resultpart[i] = _cell2region(tc)
    end

    r = zeros(Int, ndim * 2 + 1, 0)
    for i in eachindex(resultpart)
        tr = deepcopy(resultpart[i])
        tr[1:2, :] .+= limitpart[i][5] - 1
        tr[3:4, :] .+= limitpart[i][3] - 1
        tr[5:6, :] .+= limitpart[i][1] - 1
        r = hcat(r, tr)
    end

    return r
end


function fillpmlvalue!(v, mask, nx::Int, ny::Int, nz::Int, ax::Int, ay::Int, az::Int)
    segx = (1:ax, 1+ax:nx+ax, 1+nx+ax:ax+nx+ax)
    segy = (1:ay, 1+ay:ny+ay, 1+ny+ay:ay+ny+ay)
    segz = (1:nz, 1+nz:az+nz)

    # outline = [CartesianIndices((segy[2][1]:segy[2][1], segx[2]))[:];
    #            CartesianIndices((segy[2][end]:segy[2][end], segx[2]))[:];
    #            CartesianIndices((segy[2], segx[2][1]:segx[2][1]))[:];
    #            CartesianIndices((segy[2], segx[2][end]:segx[2][end]))[:]]

    for i = 1:3, j = 1:3
        if (i == 2) && (j == 2)
            continue
        end
        for jx in segx[i], jy in segy[j]
            if i == 1
                ix = segx[2][1]
            elseif i == 2
                ix = jx
            else
                ix = segx[2][end]
            end
            if j == 1
                iy = segy[2][1]
            elseif j == 2
                iy = jy
            else
                iy = segy[2][end]
            end
            for iz in segz[1]
                v[iz, jy, jx] = v[iz, iy, ix]
                mask[iz, jy, jx] = mask[iz, iy, ix]
            end
        end
    end
    # Threads.@threads for iz in segz[1]
    #     for i = 1:3, j = 1:3
    #         if (i == 2) && (j == 2)
    #             continue
    #         end
    #         for jx in segx[i], jy in segy[j]
    #             if i == 1
    #                 ix = segx[2][1]
    #             elseif i == 2
    #                 ix = jx
    #             else
    #                 ix = segx[2][end]
    #             end
    #             if j == 1
    #                 iy = segy[2][1]
    #             elseif j == 2
    #                 iy = jy
    #             else
    #                 iy = segy[2][end]
    #             end
    #             v[iz, jy, jx] = v[iz, iy, ix]
    #             mask[iz, jy, jx] = mask[iz, iy, ix]
    #         end
    #     end
    #     #
    #     # tv = v[iz, outline]
    #     # tm = mask[iz, outline]
    #     # if all(tm)
    #     #     vm = mean(tv)
    #     # else
    #     #     vm = mean(tv[.!tm])
    #     # end
    #     # for i = 1:3, j = 1:3
    #     #     if (i == 2) && (j == 2)
    #     #         continue
    #     #     end
    #     #     v[iz, segy[i], segx[j]] .= vm
    #     # end
    # end
    # for i = 1:3, j = 1:3
    #     if (i == 2) && (j == 2)
    #         continue
    #     end
    #     for jx in segx[i], jy in segy[j]
    #         if i == 1
    #             ix = segx[2][1]
    #         elseif i == 2
    #             ix = jx
    #         else
    #             ix = segx[2][end]
    #         end
    #         if j == 1
    #             iy = segy[2][1]
    #         elseif j == 2
    #             iy = jy
    #         else
    #             iy = segy[2][end]
    #         end
    #         iz = findfirst(!, mask[:, iy, ix])
    #         v[1:iz, jy, jx] .= v[iz, iy, ix]
    #     end
    # end
    # vb = mean(v[nz, :, :])
    # Threads.@threads for ix = 1:nx+2*ax
    #     v[segz[2], 1:ny+2*ay, ix] .= vb
    # end
    for iz in segz[2]
        v[iz, :, :] .= v[iz-1, :, :]
    end
    Threads.@threads for iz in segz[2]
        mask[iz, :, :] .= false
    end
    # for iz = nz:-1:1
    #     flag = true
    #     flag &= all(mask[iz, outline])
    #     for i = 1:3, j = 1:3
    #         if (i == 2) && (j == 2)
    #             continue
    #         end
    #         flag |= any(mask[iz+1, segy[i], segx[j]])
    #     end

    #     for i = 1:3, j = 1:3
    #         if (i == 2) && (j == 2)
    #             continue
    #         end
    #         mask[iz, segy[i], segx[j]] .= flag
    #     end
    # end
    return nothing
end

function appendpml(vp, vs, rho, mask, ax::Int, ay::Int, az::Int)
    (nz, ny, nx) = size(vp)
    α = zeros(nz + az, ny + 2 * ay, nx + 2 * ax)
    β = zeros(nz + az, ny + 2 * ay, nx + 2 * ax)
    ρ = zeros(nz + az, ny + 2 * ay, nx + 2 * ax)
    newmask = falses(nz + az, ny + 2 * ay, nx + 2 * ax)
    # for ix = 1:nx, iy = 1:ny, iz = 1:nz
    Threads.@threads for idx in CartesianIndices(vp)
        (iz, iy, ix) = idx.I
        α[iz, iy+ay, ix+ax] = vp[iz, iy, ix]
        β[iz, iy+ay, ix+ax] = vs[iz, iy, ix]
        ρ[iz, iy+ay, ix+ax] = rho[iz, iy, ix]
        newmask[iz, iy+ay, ix+ax] = mask[iz, iy, ix]
    end
    fillpmlvalue!(α, newmask, nx, ny, nz, ax, ay, az)
    fillpmlvalue!(β, newmask, nx, ny, nz, ax, ay, az)
    fillpmlvalue!(ρ, newmask, nx, ny, nz, ax, ay, az)
    return ((α, β, ρ), newmask)
end
