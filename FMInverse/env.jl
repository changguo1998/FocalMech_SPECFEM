PTH = Dict{String,String}()
#PTH["root"] = "/data2/guochang/chuandian_batch_rerun"
PTH["root"] = abspath(@__DIR__, "../..")
PTH["invsrc"] = abspath(@__DIR__)
PTH["invscript"] = abspath(@__DIR__, "inv_script")

PTH["raw1"] = "/data/guochang/chuandian2/dat/chuandian1"
PTH["raw2"] = "/data/guochang/chuandian2/dat/raw"
PTH["batch1"] = joinpath(PTH["root"], "dat", "unpack_batch1")
PTH["batch2"] = joinpath(PTH["root"], "dat", "unpack_batch2")
PTH["dat_filted"] = abspath(PTH["root"], "dat", "filted")
PTH["dat_selectstation"] = abspath(PTH["root"], "dat", "selectstation")
PTH["stationlist"] = joinpath(PTH["root"], "dat", "statistic", "stationlist_cb.txt")
PTH["glib"] = "/data/guochang/chuandian_green_fun/"
