julia -t 16 ../../VelocityModel/genmodel.jl \
-S 30.0/105.0 \
-D 50 \
-H 10 \
-R -25/25/-25/25 \
-N 51/51/51 \
-M CVM2.bin \
-O ./model.bin
