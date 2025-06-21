mkdir -p unpack
for f in `ls raw_data`; do
    unzip -o raw_data/"$f" ${f%.zip}.asc -d unpack
done
