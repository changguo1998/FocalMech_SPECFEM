using GLMakie

station_range = let
    buf = filter(!startswith("#"), readlines("../../dat/statistic/station_area.csv"))
    tl = map(buf) do l
        m = split(l, ' ', keepempty=false)
        (String(m[1]), parse.(Float64, m[2:end]))
    end
    tab = Matrix{Any}(undef, length(buf), 7)
    for i = eachindex(tl)
        tab[i, 1] = tl[i][1]
        tab[i, 2:7] = tl[i][2][1:6]
    end
    tab
end

downloaded_station = setdiff(readdir("../../dat/glib_hpc/"), ["index.txt", "scripts"])
didx = map(s->findfirst(==(s), station_range[:, 1]), downloaded_station)

fig = Figure(resolution=(800,1000));
ax = Axis(fig[1, 1]);
scatter!(Float64.(station_range[:, 3]), Float64.(station_range[:, 2]);
    marker=:utriangle, markersize=12, color=:gray70);
scatter!(Float64.(station_range[didx, 3]), Float64.(station_range[didx, 2]);
    marker=:utriangle, markersize=12, color=:blue);
fig
