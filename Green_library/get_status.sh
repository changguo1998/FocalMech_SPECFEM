source /etc/profile
source ~/.bashrc
if [[ -f glib_calc_status.html ]]; then
    rm glib_calc_status.html
fi
if [[ -f station_wait_to_download.txt ]]; then
    rm station_wait_to_download.txt
fi
/usr/local/bin/rayfile-c -a hefei03.nscc-hf.cn \
    -P 65012 \
    -u guochang \
    -w 9380ec0f4b63a90fbc-7701-41dc-8e42-8227f88a20d1 \
    -o download \
    -d ./ \
    -s /semsimu/chuandian2/src/glib/glib_calc_status.html && \
/usr/local/bin/rayfile-c -a hefei03.nscc-hf.cn \
    -P 65012 \
    -u guochang \
    -w 9380ec0f4b81ad476f-5ce3-41c8-b99f-da1f26cf37e2 \
    -o download \
    -d ./ \
    -s /semsimu/chuandian2/src/glib/station_wait_to_download.txt && \
/home/guochang/Softwares/julia-1.7.3/bin/julia 12.2_autodownload_local.jl