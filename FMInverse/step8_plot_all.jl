using JuliaSourceMechanism, JLD2, LinearAlgebra, Printf, CairoMakie, SeisTools

function linearscale(v::AbstractArray; ylim::AbstractVector=[-0.5, 0.5], A::Union{Real,Nothing}=nothing)
    if isnothing(A)
        A = maximum(abs, v)
    end
    y = @. (v + A) / 2 / A * (ylim[2] - ylim[1]) + ylim[1]
    return (y, A)
end

function plansubplot(xmin::Real, xmax::Real, ymin::Real, ymax::Real,
    tmin::AbstractMatrix, tmax::AbstractMatrix; xmargin::Float64=0.05, ymargin::Float64=0.01)
    trange = tmax - tmin
    tlen = maximum(trange, dims=1) |> x -> reshape(x, length(x))
    xperct = (xmax - xmin - (length(tlen) - 1) * xmargin) .* tlen ./ sum(tlen)
    yperct = (ymax - ymin - (size(tmin, 1) - 1) * ymargin) .* ones(size(tmin, 1)) ./ size(tmin, 1)
    xrange = zeros(size(tmin, 2), 2)
    yrange = zeros(size(tmin, 1), 2)
    for i = axes(tmin, 2)
        xrange[i, 2] = xmin + sum(xperct[1:i]) + (i - 1) * xmargin
        xrange[i, 1] = xrange[i, 2] - xperct[i]
    end
    for i = axes(tmin, 1)
        yrange[i, 2] = ymin + sum(yperct[1:i]) + (i - 1) * ymargin
        yrange[i, 1] = yrange[i, 2] - yperct[i]
    end
    return (xrange, reverse(yrange, dims=1), xperct[1] / tlen[1])
end

include("inv_script/plot_segment.jl")

if !(@isdefined hideaxis)
    const hideaxis = (xlabelvisible=false,
        xgridvisible=false,
        xticklabelsvisible=false,
        xticksvisible=false,
        ylabelvisible=false,
        ygridvisible=false,
        yticklabelsvisible=false,
        yticksvisible=false,
        leftspinevisible=false,
        rightspinevisible=false,
        topspinevisible=false,
        bottomspinevisible=false
    )

    const COLOR_RED = RGBf(240/255, 0/255, 86/255)
    const COLOR_BLUE = RGBf(40/255, 40/255, 238/255)
    const COLOR_WHITE = RGBf(1.0, 1.0, 1.0)
    const CMAP_RED = [COLOR_WHITE, COLOR_RED]
    const CMAP_BLUE = [COLOR_WHITE, COLOR_BLUE]
    const FONT_ = "Courier 10 Pitch"
    const BB_LINE_WIDTH = 3
    const BB_RESOLUTION = 201
    const WAVE_LINE_WIDTH = 2
    const WAVE_FONT_SIZE = 12
    const TAG_FONT_SIZE = 20
end

# if !@isdefined fname
#     fname = ARGS[1]
# end

# plot_func(fname)
# eventname = "201304200802479"
include("env.jl")

for e in readdir(PTH["dat_selectstation"])
    if isdir(joinpath(PTH["dat_selectstation"], e, "result"))
        for r in filter(endswith("jld2"), readdir(joinpath(PTH["dat_selectstation"], e, "result"), join=true))
            if isfile(replace(r, ".jld2"=>"_mt.png"))
                continue
            end
            # println(r)
            try
                plot_func(r)
            catch x
                println(x)
                # rm(r)
            end
        end
    end
end
