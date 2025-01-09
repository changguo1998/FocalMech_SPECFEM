station_list = readdir() |> x->filter(startswith("batch"), x) |> x->readlines.(x) |>
    x->sort.(x) |> x->vcat(x...)

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

for i = eachindex(station_list)
    id = findfirst(==(station_list[i]), station_range[:, 1])
    open("../slurm/$(station_list[i]).slurm", "w") do io
        println(io, """
#!/bin/bash
#SBATCH --job-name=generate_model_$(station_list[i])
#SBATCH --partition=hfacnormal01
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --output=%j.out
#SBATCH --error=%j.err

cd ../glib && \\
julia 4.1_write_setting_of_1_station.jl $(id) && \\
touch ../../dat/glib/$(station_list[i])/slicemodel.flag && \\
date""")
    end
end
