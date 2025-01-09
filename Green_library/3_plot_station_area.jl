using CairoMakie, DelimitedFiles

station_range = let
    buf = readlines("../../dat/statistic/station_area.csv")
    filter!(!startswith("#"), buf)
    t = Matrix{Any}(undef, length(buf), 7)
    for i = eachindex(buf)
        l = split(buf[i], ' '; keepempty=false)
        t[i, 1] = l[1]
        t[i, 2:end] = parse.(Float64, l[2:end])
    end
    t
end

(cvm2, _) = readdlm("../../dat/cvm2/SWChinaCVMv2.0.txt.wrst.sea_level"; header=true)
latlist = sort(unique(cvm2[:, 2]))
lonlist = sort(unique(cvm2[:, 1]))
deplist = sort(unique(cvm2[:, 3]))
depcount = zeros(Int, length(latlist), length(lonlist))
for i = axes(cvm2, 1)
    ilon = findfirst(abs.(lonlist .- cvm2[i, 1]) .< 1e-3)
    ilat = findfirst(abs.(latlist .- cvm2[i, 2]) .< 1e-3)
    idep = findfirst(abs.(deplist .- cvm2[i, 3]) .< 1e-3)
    depcount[ilat, ilon] += 1
end
filledmap = depcount .== length(deplist)


fig = Figure(resolution=(1600,2000))
ax = Axis(fig[1, 1]; autolimitaspect=1)
heatmap!(lonlist, latlist, permutedims(filledmap); colormap=[:white, :gray70])
for s = axes(station_range, 1)
    lx = Float64.(station_range[s, [6,6,7,7,6]])
    ly = Float64.(station_range[s, [4,5,5,4,4]])
    lines!(lx, ly; color=:black, linewidth=2)
end
scatter!(Float64.(station_range[:, 3]), Float64.(station_range[:, 2]);
    marker=:utriangle, markersize=32, color=:red)
save("../../dat/statistic/station_area_plot.png", fig)


mkpath("../../dat/station_area_plot/")
for s = axes(station_range, 1)
    fig2 = Figure(resolution=(1600,2000))
    ax2 = Axis(fig2[1, 1]; autolimitaspect=1.0/cosd(28.0),xticks=98:0.5:107, yticks=22:0.5:33)
    heatmap!(lonlist, latlist, permutedims(filledmap); colormap=[:white, :gray70])
    vlines!(collect(98:0.5:107), color=:gray10, width=1)
    hlines!(collect(22:0.5:33), color=:gray10, width=1)
    lx = Float64.(station_range[s, [6,6,7,7,6]])
    ly = Float64.(station_range[s, [4,5,5,4,4]])
    lines!(lx, ly; color=:blue, linewidth=4)
    scatter!(Float64[station_range[s, 3]], Float64[station_range[s, 2]];
        marker=:utriangle, markersize=32, color=:red)
    save(abspath("../../dat/station_area_plot", station_range[s, 1]*".png"), fig2)
end
