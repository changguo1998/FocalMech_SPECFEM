function plot_func(fname::String)

(env, status, result) = let
    t = JLD2.load(fname)
    (t["env"], t["status"], t["result"])
end

stationnames = filter(!startswith("info_"), collect(keys(result)))
dists = map(x -> result[x]["dist"], stationnames)
distperm = sortperm(dists)
nstation = length(stationnames)
WIDTH = 16
HDRHEIGHT = WIDTH/4
HEIGHT = nstation * WIDTH / 6 * 0.4 + HDRHEIGHT

tmin = zeros(max(nstation, 10), 6)
tmax = zeros(max(nstation, 10), 6)
for i = 1:nstation
    s = stationnames[i]
    for j = 1:6
        if mod(j - 1, 3) == 0
            cn = "E"
        elseif mod(j - 1, 3) == 1
            cn = "N"
        else
            cn = "Z"
        end
        pn = (j < 4) ? "P" : "S"
        tmin[i, j] = 0.0
        tmax[i, j] = 0.0
        if !(cn in keys(result[s]))
            continue
        end
        if !(pn in keys(result[s][cn]))
            continue
        end
        p = result[s][cn][pn]
        # if p["xcorr_weight"] > 0.0
            tmin[i, j] = min(0.0, p["xcorr_shift"], p["mt_shift"])
            tmax[i, j] = max(0.0, p["xcorr_shift"], p["mt_shift"]) + length(p["xcorr_rec"]) * p["xcorr_dt"]
        # end
    end
end

(xframe, yframe, kx) = plansubplot(0.1, 0.98, 0.0, 0.96, tmin, tmax; xmargin=0.01, ymargin=0.02)

fig = Figure(; resolution=round.(Int, (WIDTH, HEIGHT) .* 100))
ax3 = Axis(fig[2, 1:3]; hideaxis...)

lons = Float64[]
lats = Float64[]
tags = String[]
useAA = true

for n = 1:nstation
    i = distperm[n]
    s = stationnames[i]
    ista = findfirst(t -> t["network"] * "." * t["station"] == s, env["stations"])
    push!(lons, env["stations"][ista]["meta_lon"])
    push!(lats, env["stations"][ista]["meta_lat"])
    push!(tags, s)
    hl = (yframe[n, 1] + yframe[n, 2]) / 2
    text!(ax3, xframe[1, 1]-0.01, hl; text=@sprintf("%.1fkm\n%s\n%.0f°", env["stations"][ista]["base_distance"],
        stationnames[i], env["stations"][ista]["base_azimuth"]), align=(:right, :center), fontsize=TAG_FONT_SIZE)
    for j = 1:6
        p = (j > 3) ? "S" : "P"
        if mod(j - 1, 3) == 0
            c = "E"
        elseif mod(j - 1, 3) == 1
            c = "N"
        else
            c = "Z"
        end
        if !(c in keys(result[s]))
            continue
        end
        if !(p in keys(result[s][c]))
            continue
        end
        # println(s, " ", c, " ", p)
        # if result[s][c][p]["xcorr_weight"] > 0.0
            dt = result[s][c][p]["xcorr_dt"]
            # sh = result[s][c][p]["xcorr_shift"]
            sh1 = result[s][c][p]["xcorr_shift"]
            sh2 = result[s][c][p]["xcorr_shift"]
            # sh2 = result[s][c][p]["mt_shift"]
            t0 = min(0.0, sh1, sh2)
            t = range(0.0, length = length(result[s][c][p]["xcorr_rec"]), step = dt)
            (tr, A) = linearscale(result[s][c][p]["xcorr_rec"]; ylim=yframe[n, :])
            lines!(ax3, (t .- t0) .* kx .+ xframe[j, 1], tr; color=:black, linewidth=WAVE_LINE_WIDTH)
            if useAA
                (tr, _) = linearscale(result[s][c][p]["xcorr_syn"]; ylim=yframe[n, :], A=A)
            else
                (tr, _) = linearscale(result[s][c][p]["xcorr_syn"]; ylim=yframe[n, :])
            end
            lines!(ax3, (t .+ sh1 .- t0) .* kx .+ xframe[j, 1], tr; color=COLOR_BLUE, linewidth=WAVE_LINE_WIDTH)
            xc = SeisTools.DataProcess.xcorr_t(result[s][c][p]["xcorr_rec"], result[s][c][p]["xcorr_syn"],
                floor(Int, sh1/dt)-1, ceil(Int, sh1/dt)+1)
            xcv = maximum(xc.c)/norm(result[s][c][p]["xcorr_rec"])/norm(result[s][c][p]["xcorr_syn"])
            text!(ax3, xframe[j, 1], yframe[n, 1]; text=@sprintf("%.0f\n%.2fs", xcv*100, sh1), align=(:left, :baseline),
                color=COLOR_BLUE, fontsize=WAVE_FONT_SIZE)
            if useAA
                (tr, _) = linearscale(result[s][c][p]["mt_syn"]; ylim=yframe[n, :], A=A)
            else
                (tr, _) = linearscale(result[s][c][p]["mt_syn"]; ylim=yframe[n, :])
            end
            lines!(ax3, (t .+ sh2 .- t0) .* kx .+ xframe[j, 1], tr; color=COLOR_RED, linewidth=WAVE_LINE_WIDTH)
            xc = SeisTools.DataProcess.xcorr_t(result[s][c][p]["xcorr_rec"], result[s][c][p]["mt_syn"],
                floor(Int, sh2/dt)-1, ceil(Int, sh2/dt)+1)
            xcv = maximum(xc.c)/norm(result[s][c][p]["xcorr_rec"])/norm(result[s][c][p]["mt_syn"])
            text!(ax3, (t[end]+abs(sh2))*kx+xframe[j, 1], yframe[n, 1]; text=@sprintf("%.0f\n%.2fs", xcv*100, sh2),
                align=(:right, :baseline), color=COLOR_RED, fontsize=WAVE_FONT_SIZE)
            if (j <= 3) && (!iszero(result[s][c][p]["polarity_rec"]))
                pol_obs = (result[s][c][p]["polarity_rec"] > 0) ? "+" : "-"
                pol_syn = (result[s][c][p]["polarity_syn"] > 0) ? "+" : "-"
                text!(ax3, xframe[j, 1], yframe[n, 2]; text=pol_obs*"/"*pol_syn, align=(:left, :top))
            end
        # end
        if (n == nstation) && (c == "E")
            # wlen = length(result[s][c][p]["xcorr_rec"]) * result[s][c][p]["xcorr_dt"]
            wlen = (xframe[j+1, 1] - xframe[j, 1] - 0.01) / kx
            basepow = 10^floor(log10(wlen))
            rulerlen = floor(wlen / basepow) * basepow
            kruler = rulerlen * kx
            lines!(ax3, [xframe[j, 1], xframe[j, 1] + kruler], [1.0, 1.0].*yframe[n, 1]; color=:black)
            lines!(ax3, [1.0, 1.0].*xframe[j, 1], [0.0, -0.1].*(yframe[n, 2]-yframe[n, 1]).+yframe[n, 1]; color=:black)
            lines!(ax3, [1.0, 1.0].*(xframe[j, 1] + kruler),
                [0.0, -0.1].*(yframe[n, 2]-yframe[n, 1]).+yframe[n, 1]; color=:black)
            text!(ax3, xframe[j, 1] + kruler / 2, yframe[n, 1]; text = @sprintf("%gs", rulerlen),
                align=(:center, :top))
        end
    end
end
text!(ax3, (xframe[:, 1] .+ xframe[:, 2]) ./ 2, ones(6);
    text = ["PE", "PN", "PZ", "SE", "SN", "SZ"], align=(:center, :top), fontsize=TAG_FONT_SIZE)
xlims!(ax3, 0.0, 1.0)
ylims!(ax3, 0.0, 1.0)

ax2 = Axis(fig[1, 3], autolimitaspect=cosd(env["event"]["latitude"]))
scatter!(ax2, lons, lats, marker=:utriangle, color=:blue, markersize=25)
scatter!(ax2, [env["event"]["longitude"]], [env["event"]["latitude"]], marker=:star5, color=:red, markersize=36)
text!(ax2, lons, lats; text=tags, align=(:center, :top))


ax1 = Axis(fig[1, 1]; autolimitaspect=1, limits=(-1.1, 1.1, -3.1, 1.1), hideaxis...)
mt = SeisTools.Source.MomentTensor(result["info_mech"][1], result["info_mech"][2], result["info_mech"][3])
mt_plane = SeisTools.Source.focalmechanism(mt)
x0 = range(; start=-1.0, stop=1.0, length=BB_RESOLUTION)
heatmap!(ax1, x0, x0,
    SeisTools.Source.beachball_bitmap(mt; resolution=(BB_RESOLUTION, BB_RESOLUTION)) |>
    permutedims .|> sign;
    colormap=CMAP_BLUE)
l = SeisTools.Source.beachball_sdrline(mt)
lines!(ax1, map(x -> (x[2], x[1]), l.l1); color=:black, linewidth=2)
lines!(ax1, map(x -> (x[2], x[1]), l.l2); color=:black, linewidth=2)
lines!(ax1, map(x -> (x[2], x[1]), l.edge); color=:black, linewidth=4)
text!(ax1, 0.0, -1.0, text=@sprintf("longitude:%g\nlatutude:%g\ndepth:%g km\n   Mw: %.1f\nplane1: %3d/%2d/%3d\nplane2: %3d/%2d/%3d",
    env["event"]["longitude"], env["event"]["latitude"], env["algorithm"]["searchdepth"], result["info_mag"],
    mt_plane.plane1..., mt_plane.plane2...),
    align=(:center, :top), fontsize=26, font=FONT_)


ax4 = Axis(fig[1, 2]; autolimitaspect=1, limits=(-1.1, 1.1, -3.1, 1.1), hideaxis...)
_tmt = SeisTools.Source.MomentTensor(result["info_mt"])
tmt = SeisTools.Source.MomentTensor(result["info_mt"]./SeisTools.Source.Fnorm(_tmt).*sqrt(2))
mt_dcp = SeisTools.Source.decompose(tmt)
tmt_plane = SeisTools.Source.focalmechanism(mt_dcp.dc)
heatmap!(ax4, x0, x0,
    SeisTools.Source.beachball_bitmap(tmt; resolution=(BB_RESOLUTION, BB_RESOLUTION)) |>
    permutedims .|> sign;
    colormap=CMAP_RED,
    colorrange=(-1.0, 1.0))
l2 = SeisTools.Source.beachball_sdrline(tmt)
lines!(ax4, map(x -> (x[2], x[1]), l2.l1); color=:black, linewidth=2)
lines!(ax4, map(x -> (x[2], x[1]), l2.l2); color=:black, linewidth=2)
lines!(ax4, map(x -> (x[2], x[1]), l2.edge); color=:black, linewidth=4)
lines!(ax4, map(x -> (x[2], x[1]), l.l1); color=:black, linewidth=2, linestyle=:dash)
lines!(ax4, map(x -> (x[2], x[1]), l.l2); color=:black, linewidth=2, linestyle=:dash)
text!(ax4, 0.0, -1.0, text=@sprintf("plane1: %3d/%2d/%3d\nplane2: %3d/%2d/%3d\nDC   %5.1f%%\nISO  %5.1f%%\nCLVD %5.1f%%",
        tmt_plane.plane1..., tmt_plane.plane2..., 50*SeisTools.Source.Fnorm(mt_dcp.dc)^2,
        50*SeisTools.Source.Fnorm(mt_dcp.iso)^2, 50*SeisTools.Source.Fnorm(mt_dcp.clvd)^2),
        align=(:center, :top), font=FONT_, fontsize=26)

rowsize!(fig.layout, 1, Fixed(round(Int, HDRHEIGHT*100)))
colsize!(fig.layout, 1, Aspect(1, 1.0))
colsize!(fig.layout, 2, Aspect(1, 1.0))
save(replace(fname, ".jld2" => "_mt.png"), fig; px_per_unit=4)
# save(replace(fname, ".jld2" => "_mt.pdf"), fig; px_per_unit=4)

fig = Figure()
ax = Axis(fig[1, 1])
mvsd_d = result["info_misvsdep"][1]
mvsd_m = result["info_misvsdep"][2]
scatterlines!(mvsd_d[sortperm(mvsd_d)], mvsd_m[sortperm(mvsd_d)])
save(replace(fname, ".jld2" => "_mis_vs_dep.png"), fig)
# save(replace(fname, ".jld2" => "_mis_vs_dep.pdf"), fig)

end