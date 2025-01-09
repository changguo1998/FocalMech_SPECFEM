using ArgumentProcessor, TOML, Printf

include("../include/SEMshell.jl")
include("../include/io.jl")

addopt!("settingfile"; abbr="S", default=" " * joinpath(pwd(), "setting.toml"),
        fmt=" %s", help="path to setting file")
addopt!("expbit"; abbr="Be", fmt=" %d", default=" 3", help="nbit of exponent")
addopt!("sigbit"; abbr="Bs", fmt=" %d", default=" 7", help="nbit of significand")
addopt!("nthread"; abbr="T", fmt=" %d", default=" 16", help="nbit of significand")

function main()
    input = ArgumentProcessor.parse(ARGS)
    compressExe = abspath(@__DIR__, "../semexec/bin/compress")
    # - constants
    shotsetting = SEMshell.loadshotsetting(abspath(input.settingfile))
    # - compress glib files
    for ilib = eachindex(shotsetting[4].receiver)
        println("compress ", ilib, " of ", length(shotsetting[4].receiver))
        lrec = shotsetting[4].receiver[ilib]
        println("  id: ", lrec.id)
        # - file path
        inputpath = abspath(shotsetting[1].wkdir, "..", @sprintf("glib_tmp1_%d.bin", lrec.id))
        glibpath  = abspath(shotsetting[1].wkdir, "..", @sprintf("glib_%d.bin", lrec.id))
        run(pipeline(Cmd(
            Cmd([
                compressExe, inputpath, glibpath, string(input.expbit), 
                string(input.sigbit), string(input.nthread)
                ]);
            dir=abspath(shotsetting[1].wkdir, "..")
        ), devnull))
    end
end

main()
