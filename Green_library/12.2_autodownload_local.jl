fs = readlines("station_wait_to_download.txt")

if isempty(fs)
    exit(0)
end

if isfile("autodownload.flag")
    exit(0)
end

touch("autodownload.flag")

localglibdir = abspath("../../dat/glib_hpc/")
remoteglibdir= "/semsimu/chuandian2/dat/glib"
gstation = fs[1]
cmd_update_downloading_flag = Cmd([
    "/usr/local/bin/rayfile-c",
    "-a", "hefei03.nscc-hf.cn",
    "-P", "65012",
    "-u", "guochang",
    "-w", "b1ed11352feac6143e-f7ea-446c-9303-bca6b3bd9de1",
    "-o", "upload",
    "-d", joinpath(remoteglibdir, gstation),
    "-s", "downloading.flag"])

cmd_download_data = Cmd([
    "/usr/local/bin/rayfile-c",
    "-a", "hefei03.nscc-hf.cn",
    "-P", "65012",
    "-u", "guochang",
    "-w", "9380ec0f4b81ad476f-5ce3-41c8-b99f-da1f26cf37e2",
    "-o", "download",
    "-d", localglibdir,
    "-s", joinpath(remoteglibdir, gstation)])

cmd_update_download_flag = Cmd([
    "/usr/local/bin/rayfile-c",
    "-a", "hefei03.nscc-hf.cn",
    "-P", "65012",
    "-u", "guochang",
    "-w", "b1ed11352feac6143e-f7ea-446c-9303-bca6b3bd9de1",
    "-o", "upload",
    "-d", joinpath(remoteglibdir, gstation),
    "-s", "download.flag"])

run(Cmd(cmd_update_downloading_flag; dir=@__DIR__))
run(Cmd(cmd_download_data; dir=@__DIR__))
run(Cmd(cmd_update_download_flag; dir=@__DIR__))

rm("autodownload.flag")