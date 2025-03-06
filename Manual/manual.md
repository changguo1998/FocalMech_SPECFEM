# VelocityModel

Scripts in this directory is used to preprocess 3D velocity model, which will be
used in the calculation of 3D Green's Function.

## Data Downloading

### Elevation data

Users may want to download elevation and 3D velocity data by themselves.

The digital elevation model (DEM) can be found in [CGIAR-SRTM webpage](https://srtm.csi.cgiar.org/srtmdata/)

### 3D models:

- [SWChinaCVM1.0](https://doi.org/10.1785/0220200318): Ying Liu, Huajian Yao, Haijiang Zhang, Hongjian Fang; The Community Velocity Model V.1.0 of Southwest China, Constructed from Joint Body‐ and Surface‐Wave Travel‐Time Tomography. Seismological Research Letters 2021; 92 (5): 2972–2987.
- [SWChinaCVM2.0](https://doi.org/10.1007/s11430-022-1161-7): Liu, Y., Yu, Z., Zhang, Z. et al. The high-resolution community velocity model V2.0 of southwest China, constructed by joint body and surface wave tomography of data recorded at temporary dense arrays. Sci. China Earth Sci. 66, 2368–2385 (2023).

## Data format

### 1. binary velocity model

| type                 | variable name | meaning                |
| :------------------- | :-----------: | :--------------------- |
| int32                |      nx       | points along latitude  |
| int32                |      ny       | points along longitude |
| int32                |      nz       | points along depth     |
| float32 * `nx`       |     lats      | latitudes              |
| float32 * `ny`       |     lons      | longitudes             |
| float32 * `nz`       |     deps      | depth                  |
| float32 * `nx*ny*nz` |      vp       | P velocity             |
| float32 * `nx*ny*nz` |      vs       | S velocity             |

### 2. binary topography data

| type                  | variable name | meaning                |
| :-------------------- | :-----------: | :--------------------- |
| int32                 |     nlon      | points along longitude |
| int32                 |     nlat      | points along latitude  |
| float32 * `nlon`      |     lons      | longitudes             |
| float32 * `nlat`      |     lats      | latitudes              |
| float32 * `nlon*nlat` |     topo      | topography data        |

## Scripts

### `convertmodel_txt2bin.jl`

This script preprocesses `SWChinaCVM-2.0` model. Users need to download SWChinaCVM2.0 model and unpack it in current directory,
then run this script to get prepossed binary file.

### `converttopodata.jl`

This script preprocess SRTM DEM data. Uses need to download 5x5 Esri ascii DEM data, and unpack all to `unpack` directory,
then run this script to get combined binary file. Users may need to modify the index according to
downloaded block.

### `genmodel.jl`

This scripts is used to cut 3d model with topography from total model. Users can run
`genmodel.jl -h` to see explainations for parameters
