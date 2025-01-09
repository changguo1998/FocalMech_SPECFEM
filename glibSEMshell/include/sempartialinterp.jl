const Int_ = Int32
const Float_ = Float32

if !@isdefined(MAX_NETWORK_LEN)
    const HEADER_LEN      = 512
    const MAX_NETWORK_LEN = 8
    const MAX_STATION_LEN = 32
end

macro _sempartial_checksize(var, vsize)
    return :(if $(esc(:(size($var) != $vsize)))
                 error("SizeError: size of " * $(string(var)) * " must be " * $(esc(:(string($vsize)))))
             end)
end

function _readarray(io::IO, T::Type, s::Tuple)
    x = zeros(T, s)
    read!(io, x)
    return x
end

function _readarray(io::IO, T::Type, n::Integer)
    return _readarray(io, T, (n,))
end

function readglob(io::IO)
    c = _readarray(io, Int_, 8)
    (NDIM, NGLOB, NSPEC, NT, NGLLX, NGLLY, NGLLZ, NREC) = c
    if (NGLLX != NGLLY) || (NGLLY != NGLLZ)
        error("NGLL? not equal")
    end
    NGLL = NGLLX
    spec2glob = _readarray(io, Int_, (NGLL, NGLL, NGLL, NSPEC))
    rec2spec = _readarray(io, Int_, NREC)
    l = zeros(Float_, NDIM, NDIM, NGLL, NGLL, NGLL, NSPEC)
    begin
        read!(io, @view(l[1, 1, :, :, :, :])) # ξ₁
        read!(io, @view(l[2, 1, :, :, :, :])) # ξ₂
        read!(io, @view(l[3, 1, :, :, :, :])) # ξ₃
        read!(io, @view(l[1, 2, :, :, :, :])) # η₁
        read!(io, @view(l[2, 2, :, :, :, :])) # η₂
        read!(io, @view(l[3, 2, :, :, :, :])) # η₃
        read!(io, @view(l[1, 3, :, :, :, :])) # γ₁
        read!(io, @view(l[2, 3, :, :, :, :])) # γ₂
        read!(io, @view(l[3, 3, :, :, :, :])) # γ₃
    end
    tnet    = _readarray(io, UInt8, (MAX_NETWORK_LEN, NREC))
    network = map(i -> String(tnet[:, i]), 1:NREC)
    tsta    = _readarray(io, UInt8, (MAX_STATION_LEN, NREC))
    station = map(i -> String(tsta[:, i]), 1:NREC)
    hpxx    = _readarray(io, Float_, (NGLL, NGLL))
    hpyy    = _readarray(io, Float_, (NGLL, NGLL))
    hpzz    = _readarray(io, Float_, (NGLL, NGLL))
    hp      = zeros(Float_, NDIM, NGLL, NGLL)
    for i = 1:NGLLX, j = 1:NGLLX
        hp[1, i, j] = hpxx[j, i]
        hp[2, i, j] = hpyy[j, i]
        hp[3, i, j] = hpzz[j, i]
    end
    νr = _readarray(io, Float64, (NDIM, NDIM, NREC))
    return (NDIM, NGLOB, NSPEC, NT, NGLL, NREC, spec2glob, rec2spec, l, network, station, hp, νr)
end

function readlocalhead(io::IO)
    c = _readarray(io, Int_, 5)
    (NDIM, NGLL, NT, LREC, LGLOB) = c
    lglob2glob = _readarray(io, Int_, LGLOB)
    lrec2rec = _readarray(io, Int_, LREC)
    ξr = _readarray(io, Float_, (NGLL, LREC))
    ηr = _readarray(io, Float_, (NGLL, LREC))
    γr = _readarray(io, Float_, (NGLL, LREC))
    return (NDIM, NT, LREC, LGLOB, lglob2glob, lrec2rec, ξr, ηr, γr)
end

function readlocalfield(io::IO)
    (NDIM, NT, LGLOB) = _readarray(io, Int_, 3)
    v = _readarray(io, Float_, (NDIM, LGLOB, NT))
    return (NDIM, NT, LGLOB, v)
end

function prepare_interpfield(NDIM::Integer, LREC::Integer, LGLOB::Integer, lglob2glob::AbstractVector{Int_})
    glob2lglob = Dict{Int_,Int}()
    for i = 1:LGLOB
        glob2lglob[lglob2glob[i]] = i
    end
    tmp = zeros(Float_, NDIM, LREC)
    return (tmp, glob2lglob)
end

function interpolatefield!(interped::AbstractMatrix{Float_}, field::AbstractMatrix{Float_},
                           NDIM::Int_, NSPEC::Int_, NGLL::Int_, NREC::Int_, rec2spec::AbstractVector{Int_},
                           spec2glob::AbstractArray{Int_,4}, νr::AbstractArray{Float64,3},
                           LREC::Int_, LGLOB::Int_, lrec2rec::AbstractVector{Int_},
                           ξr::AbstractMatrix{Float_}, ηr::AbstractMatrix{Float_}, γr::AbstractMatrix{Float_},
                           glob2lglob::Dict{Int_,Int}, tmp::AbstractMatrix{Float_})
    # * * * * * * * * * *
    # *   check size    *
    # * * * * * * * * * *
    @_sempartial_checksize interped (NDIM, LREC)
    @_sempartial_checksize field (NDIM, LGLOB)
    @_sempartial_checksize rec2spec (NREC,)
    @_sempartial_checksize spec2glob (NGLL, NGLL, NGLL, NSPEC)
    @_sempartial_checksize νr (NDIM, NDIM, NREC)
    @_sempartial_checksize lrec2rec (LREC,)
    # @_sempartial_checksize lglob2glob (LGLOB,)
    @_sempartial_checksize ξr (NGLL, LREC)
    @_sempartial_checksize ηr (NGLL, LREC)
    @_sempartial_checksize γr (NGLL, LREC)
    @_sempartial_checksize tmp (NDIM, LREC)
    #=
        glob2lglob = Dict{Int_,Int}()
        for i = 1:LGLOB
            glob2lglob[lglob2glob[i]] = i
        end
    =#
    # * * * * * * * * * *
    # *   Interpolate   *
    # * * * * * * * * * *
    # tmp = zeros(Float_, NDIM, LREC)
    Threads.@threads for irl = 1:LREC
        ir = lrec2rec[irl]
        ispec = rec2spec[ir]
        for i = 1:NDIM
            tmp[i, irl] = 0.0
            for iγ = 1:NGLL, iη = 1:NGLL, iξ = 1:NGLL
                iglob = spec2glob[iξ, iη, iγ, ispec]
                lglob = glob2lglob[iglob]
                tmp[i, irl] += field[i, lglob] * ξr[iξ, irl] * ηr[iη, irl] * γr[iγ, irl]
            end
        end
        for i = 1:NDIM
            interped[i, irl] = 0.0
            for j = 1:NDIM
                interped[i, irl] += νr[i, j, ir] * tmp[j, irl]
            end
        end
    end
    return nothing
end

function prepare_∇rec(NDIM::Integer, NGLL::Integer, LREC::Integer, LGLOB::Integer, lglob2glob::AbstractVector{Int_})
    tmp  = zeros(Float_, NDIM, LREC)
    tmp1 = zeros(Float_, NGLL, NGLL, NGLL, LREC)
    tmp2 = zeros(Float_, NDIM, NGLL, NGLL, NGLL, LREC)
    tmp3 = zeros(Float_, NDIM, NGLL, NGLL, NGLL, LREC)

    glob2lglob = Dict{Int_,Int}()
    for i = 1:LGLOB
        glob2lglob[lglob2glob[i]] = i
    end
    return (tmp, tmp1, tmp2, tmp3, glob2lglob)
end

function ∇andinterp_looprec!(interped::AbstractArray{Float_,2}, f::AbstractVector{Float_},
                             NDIM::Int_, NSPEC::Int_, NGLL::Int_, NREC::Int_, rec2spec::AbstractVector{Int_},
                             spec2glob::AbstractArray{Int_,4}, νr::AbstractArray{Float64,3},
                             hp::AbstractArray{Float_,3}, ℓ::AbstractArray{Float_,6},
                             LREC::Int_, LGLOB::Int_, lrec2rec::AbstractVector{Int_},
                             ξr::AbstractMatrix{Float_}, ηr::AbstractMatrix{Float_}, γr::AbstractMatrix{Float_},
                             glob2lglob::Dict{Int_,Int}, tmp::AbstractMatrix{Float_}, tmp1::AbstractArray{Float_,4},
                             tmp2::AbstractArray{Float_,5}, tmp3::AbstractArray{Float_,5})
    # * * * * * * * * * *
    # *   check size    *
    # * * * * * * * * * *
    @_sempartial_checksize interped (NDIM, LREC)
    @_sempartial_checksize f (LGLOB,)
    @_sempartial_checksize rec2spec (NREC,)
    @_sempartial_checksize spec2glob (NGLL, NGLL, NGLL, NSPEC)
    # @_sempartial_checksize spec2irspec (NSPEC,)
    @_sempartial_checksize νr (NDIM, NDIM, NREC)
    @_sempartial_checksize hp (NDIM, NGLL, NGLL)
    @_sempartial_checksize ℓ (NDIM, NDIM, NGLL, NGLL, NGLL, NSPECIR)
    @_sempartial_checksize lrec2rec (LREC,)
    # @_sempartial_checksize lglob2glob (LGLOB,)
    @_sempartial_checksize ξr (NGLL, LREC)
    @_sempartial_checksize ηr (NGLL, LREC)
    @_sempartial_checksize γr (NGLL, LREC)
    @_sempartial_checksize tmp (NDIM, LREC)
    @_sempartial_checksize tmp1 (NGLL, NGLL, NGLL, LREC)
    @_sempartial_checksize tmp2 (NDIM, NGLL, NGLL, NGLL, LREC)
    @_sempartial_checksize tmp3 (NDIM, NGLL, NGLL, NGLL, LREC)
    # * * * * * * * * * *
    # *   build index   *
    # * * * * * * * * * *
    #=
    tmp  = zeros(Float_, NDIM, LREC)
    tmp1 = zeros(Float_, NGLL, NGLL, NGLL, LREC)
    tmp2 = zeros(Float_, NDIM, NGLL, NGLL, NGLL, LREC)
    tmp3 = zeros(Float_, NDIM, NGLL, NGLL, NGLL, LREC)

    glob2lglob = Dict{Int_,Int}()
    for i = 1:LGLOB
        glob2lglob[lglob2glob[i]] = i
    end
    =#
    # * * * * * * * * * *
    # *    calculate    *
    # * * * * * * * * * *
    Threads.@threads for irl = 1:LREC
        ir = lrec2rec[irl]
        ispec = rec2spec[ir]
        # ispecir = spec2irspec[ispec]
        # * * * * * * * * * *
        # *    partial      *
        # * * * * * * * * * *
        for k = 1:NGLL, j = 1:NGLL, i = 1:NGLL
            iglob = spec2glob[i, j, k, ispec]
            lglob = glob2lglob[iglob]
            tmp1[i, j, k, irl] = f[lglob]
        end
        for k = 1:NGLL, j = 1:NGLL, i = 1:NGLL
            tmp2[:, i, j, k, irl] .= 0.0
            for l = 1:NGLL
                tmp2[1, i, j, k, irl] += tmp1[l, j, k, irl] * hp[1, l, i]
                tmp2[2, i, j, k, irl] += tmp1[i, l, k, irl] * hp[2, l, j]
                tmp2[3, i, j, k, irl] += tmp1[i, j, l, irl] * hp[3, l, k]
            end
        end
        tmp3[:, :, :, :, irl] .= 0.0
        for k = 1:NGLL, j = 1:NGLL, i = 1:NGLL, idim = 1:NDIM, jdim = 1:NDIM
            tmp3[idim, i, j, k, irl] += ℓ[idim, jdim, i, j, k, ispec] * tmp2[jdim, i, j, k, irl]
        end
        # * * * * * * * * * *
        # *   interpolate   *
        # * * * * * * * * * *
        for i = 1:NDIM
            tmp[i, irl] = 0.0
            for iγ = 1:NGLL, iη = 1:NGLL, iξ = 1:NGLL
                tmp[i, irl] += tmp3[i, iξ, iη, iγ, irl] * ξr[iξ, irl] * ηr[iη, irl] * γr[iγ, irl]
            end
        end
        # rotate
        for i = 1:NDIM
            interped[i, irl] = 0.0
            for j = 1:NDIM
                interped[i, irl] += νr[i, j, ir] * tmp[j, irl]
            end
        end
    end
    return nothing
end

function prepare_∇spec(NDIM::Integer, NSPEC::Integer, NGLL::Integer, LREC::Integer, lrec2rec::Vector{Int_},
                       LGLOB::Integer,
                       lglob2glob::AbstractVector{Int_}, NREC::Integer, rec2spec::Vector{Int_})
    @_sempartial_checksize rec2spec (NREC,)
    @_sempartial_checksize lrec2rec (LREC,)
    @_sempartial_checksize lglob2glob (LGLOB,)
    glob2lglob = Dict{Int_,Int}()
    for i = 1:LGLOB
        glob2lglob[lglob2glob[i]] = i
    end

    specused = falses(NSPEC)
    for irl = 1:LREC
        ir = lrec2rec[irl]
        ispec = rec2spec[ir]
        specused[ispec] = true
    end
    LSPEC = count(specused)
    lspec2spec = findall(specused)
    spec2lspec = Dict{Int,Int}()
    for lspec = 1:LSPEC
        spec2lspec[lspec2spec[lspec]] = lspec
    end

    tmp  = zeros(Float_, NDIM, LREC)
    tmp1 = zeros(Float_, NGLL, NGLL, NGLL, LSPEC)
    tmp2 = zeros(Float_, NDIM, NGLL, NGLL, NGLL, LSPEC)
    tmp3 = zeros(Float_, NDIM, NGLL, NGLL, NGLL, LSPEC)
    return (tmp, tmp1, tmp2, tmp3, LSPEC, glob2lglob, lspec2spec, spec2lspec)
end

function ∇andinterp_loopspec!(interped::AbstractArray{Float_,2}, f::AbstractVector{Float_},
                              NDIM::Int_, NSPEC::Int_, NGLL::Int_, NREC::Int_, rec2spec::AbstractVector{Int_},
                              spec2glob::AbstractArray{Int_,4}, νr::AbstractArray{Float64,3},
                              hp::AbstractArray{Float_,3}, ℓ::AbstractArray{Float_,6},
                              LREC::Int_, LGLOB::Int_, lrec2rec::AbstractVector{Int_},
                              ξr::AbstractMatrix{Float_}, ηr::AbstractMatrix{Float_}, γr::AbstractMatrix{Float_},
                              glob2lglob::Dict{Int_,Int}, spec2lspec::Dict{Int,Int}, LSPEC::Integer,
                              lspec2spec::Vector{Int}, tmp::AbstractMatrix{Float_}, tmp1::AbstractArray{Float_,4},
                              tmp2::AbstractArray{Float_,5}, tmp3::AbstractArray{Float_,5})
    # * * * * * * * * * *
    # *   check size    *
    # * * * * * * * * * *
    @_sempartial_checksize interped (NDIM, LREC)
    @_sempartial_checksize f (LGLOB,)
    @_sempartial_checksize rec2spec (NREC,)
    @_sempartial_checksize spec2glob (NGLL, NGLL, NGLL, NSPEC)
    @_sempartial_checksize νr (NDIM, NDIM, NREC)
    @_sempartial_checksize hp (NDIM, NGLL, NGLL)
    @_sempartial_checksize ℓ (NDIM, NDIM, NGLL, NGLL, NGLL, NSPEC)
    @_sempartial_checksize lrec2rec (LREC,)
    @_sempartial_checksize ξr (NGLL, LREC)
    @_sempartial_checksize ηr (NGLL, LREC)
    @_sempartial_checksize γr (NGLL, LREC)
    @_sempartial_checksize tmp (NDIM, LREC)
    @_sempartial_checksize tmp1 (NGLL, NGLL, NGLL, LSPEC)
    @_sempartial_checksize tmp2 (NDIM, NGLL, NGLL, NGLL, LSPEC)
    @_sempartial_checksize tmp3 (NDIM, NGLL, NGLL, NGLL, LSPEC)
    # * * * * * * * * * *
    # *   build index   *
    # * * * * * * * * * *
    #=
        glob2lglob = Dict{Int_,Int}()
        for i = 1:LGLOB
            glob2lglob[lglob2glob[i]] = i
        end

        specused = falses(NSPEC)
        for irl = 1:LREC
            ir = lrec2rec[irl]
            ispec = rec2spec[ir]
            specused[ispec] = true
        end
        LSPEC = count(specused)
        lspec2spec = findall(specused)
        spec2lspec = Dict{Int,Int}()
        for lspec = 1:LSPEC
            spec2lspec[lspec2spec[lspec]] = lspec
        end

        tmp  = zeros(Float_, NDIM, LREC)
        tmp1 = zeros(Float_, NGLL, NGLL, NGLL, LSPEC)
        tmp2 = zeros(Float_, NDIM, NGLL, NGLL, NGLL, LSPEC)
        tmp3 = zeros(Float_, NDIM, NGLL, NGLL, NGLL, LSPEC)
    =#
    # * * * * * * * * * *
    # *    partial      *
    # * * * * * * * * * *
    Threads.@threads for lspec = 1:LSPEC
        ispec = lspec2spec[lspec]
        for k = 1:NGLL, j = 1:NGLL, i = 1:NGLL
            iglob = spec2glob[i, j, k, ispec]
            lglob = glob2lglob[iglob]
            tmp1[i, j, k, lspec] = f[lglob]
        end
        tmp2[:, :, :, :, lspec] .= 0.0
        for k = 1:NGLL, j = 1:NGLL, i = 1:NGLL, l = 1:NGLL
            tmp2[1, i, j, k, lspec] += tmp1[l, j, k, lspec] * hp[1, l, i]
            tmp2[2, i, j, k, lspec] += tmp1[i, l, k, lspec] * hp[2, l, j]
            tmp2[3, i, j, k, lspec] += tmp1[i, j, l, lspec] * hp[3, l, k]
        end
        tmp3[:, :, :, :, lspec] .= 0.0
        for k = 1:NGLL, j = 1:NGLL, i = 1:NGLL, idim = 1:NDIM, jdim = 1:NDIM
            tmp3[idim, i, j, k, lspec] += ℓ[idim, jdim, i, j, k, ispec] * tmp2[jdim, i, j, k, lspec]
        end
    end

    # * * * * * * * * * *
    # *   interpolate   *
    # * * * * * * * * * *
    Threads.@threads for irl = 1:LREC
        ir = lrec2rec[irl]
        ispec = rec2spec[ir]
        lspec = spec2lspec[ispec]
        tmp[:, irl] .= 0.0
        for i = 1:NDIM, iγ = 1:NGLL, iη = 1:NGLL, iξ = 1:NGLL
            tmp[i, irl] += tmp3[i, iξ, iη, iγ, lspec] * ξr[iξ, irl] * ηr[iη, irl] * γr[iγ, irl]
        end
        # rotate
        interped[:, irl] .= 0.0
        for i = 1:NDIM, j = 1:NDIM
            interped[i, irl] += νr[i, j, ir] * tmp[j, irl]
        end
    end

    return nothing
end
