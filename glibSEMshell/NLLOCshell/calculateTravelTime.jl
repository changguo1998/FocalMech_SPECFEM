using ArgumentProcessor, TOML

include("NLLOCshell.jl")

SGlobal = NLLOCshell.NLLocSettingGlobal
SVG = NLLOCshell.NLLocSettingVel2Grid3D
SGT = NLLOCshell.NLLocSettingGrid2Time
SSRC = NLLOCshell.NLLocSettingSourceLocation
ReadTT = NLLOCshell.readgrid

addopt!("model"; abbr="M", fmt=" %s", required=true, help="Path to fd model")
addopt!("wkdir"; abbr="C", fmt=" %s", default=" " * pwd(), help="Path to store temporary file")
addopt!("ttlib"; abbr="T", fmt=" %s", default=" " * joinpath(pwd(), "default.tt"),
        help="Path to output traveltime field file")
addopt!("setting"; abbr="S", fmt=" %s", help="Path to setting.toml file")
addopt!("refdx"; abbr="H", fmt=" %g", default=" 0.0", help="reference dx")
input = ArgumentProcessor.parse(ARGS)
inpath = input.model
wkdir = input.wkdir
target = input.ttlib
setting = TOML.parsefile(abspath(input.setting))

# inpath = abspath("changning.bin")
# wkdir = abspath("var")
# target = abspath("TEST.tt")
# setting = TOML.parsefile("setting.toml")
@info "Load model"
inmodel = NLLOCshell.Model(inpath)
if !((inmodel.dx == inmodel.dy) && (inmodel.dy == inmodel.dz))
    @info "Resample model"
    # inmodel = NLLOCshell.resample_nearest(inmodel, min(inmodel.dx, inmodel.dy, inmodel.dz)*2)
    inmodel = NLLOCshell.resample_nearest(inmodel, max(input.refdx, min(inmodel.dx, inmodel.dy, inmodel.dz)))
end
GC.gc()
mkpath(wkdir)
execdir = abspath(@__DIR__, "nllocexec")
for phase in ('P', 'S')
    @info "Calc $(phase) traveltime"
    phasewkdir = abspath(wkdir, String([phase]))
    velfilename = joinpath(phasewkdir, "velocity." * phase)
    mkpath(phasewkdir)
    if phase == 'P'
        open(velfilename, "w") do io
            NLLOCshell.writemodel(io, inmodel, 1)
        end
    else
        open(velfilename, "w") do io
            NLLOCshell.writemodel(io, inmodel, 2)
        end
    end
    glob = SGlobal(; trans=NLLOCshell.TransformType(; name="NONE"))
    v2g = SVG(; inputpath=velfilename,
              inputfiletype="FDTOMO",
              output=joinpath(phasewkdir, "velgrid.bin"),
              wavetypes=[String([phase])],
              nxi=size(inmodel.index, 2),
              nyi=size(inmodel.index, 3),
              nzi=size(inmodel.index, 1),
              oxi=0.0,
              oyi=0.0,
              ozi=inmodel.zt,
              dxi=inmodel.dy,
              dyi=inmodel.dx,
              dzi=inmodel.dz,
              vmin=0.0,
              vmax=10.0,
              nx=size(inmodel.index, 2),
              ny=size(inmodel.index, 3),
              nz=size(inmodel.index, 1),
              ox=0.0,
              oy=0.0,
              oz=inmodel.zt,
              dx=inmodel.dy,
              dy=inmodel.dx,
              dz=inmodel.dz,
              type="SLOW_LEN")
    g2t = SGT(; gridfile=joinpath(phasewkdir, "velgrid.bin"),
              ttfile=joinpath(phasewkdir, "ttfield.bin"),
              wavetypes=[String([phase])],
              gridmode="GRID3D",
              anglemode="ANGLES_NO",
              source=[SSRC(; label="SRC",
                           coordinate="XYZ",
                           paraf=[setting["srcloc"][2], setting["srcloc"][1]],
                           z=setting["srcloc"][3])])
    open(joinpath(phasewkdir, "nlloc.in"), "w") do io
        NLLOCshell.printsetting(io, [glob, v2g, g2t])
    end
    @info "  Run Vel2Grid3D"
    run(Cmd(Cmd([joinpath(execdir, "Vel2Grid3D"), "nlloc.in"]); dir=phasewkdir))
    @info "  Run Grid2Time"
    run(Cmd(Cmd([joinpath(execdir, "Grid2Time"), "nlloc.in"]); dir=phasewkdir))
end

GC.gc()

TP = ReadTT(joinpath(wkdir, "P", "ttfield.bin.P.SRC.time"))
TS = ReadTT(joinpath(wkdir, "S", "ttfield.bin.S.SRC.time"))
tt = zeros(Float32, 2, size(TP.data)...)
tt[1, :, :, :] .= TP.data
tt[2, :, :, :] .= TS.data

let
    nx = size(TP.data, 3)
    ny = size(TP.data, 2)
    nz = size(TP.data, 1)
    dx = TP.delta[1]
    dy = TP.delta[2]
    dz = TP.delta[3]
    ox = TP.orig[1]
    oy = TP.orig[2]
    oz = TP.orig[3]
    open(target, "w") do io
        write(io, Int32.([nx, ny, nz]))
        write(io, Float32.([dx, dy, dz, ox-setting["srcloc"][1], oy-setting["srcloc"][2], oz]))
        write(io, tt)
    end
end
