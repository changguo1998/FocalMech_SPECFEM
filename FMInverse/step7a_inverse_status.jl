include("env.jl")

checksym = "&#x2713;"

if !isfile("batchall_inverse_log.txt")
    errorevent = String[]
else
    errorevent = readlines("batchall_inverse_log.txt")
end

events = readdir(PTH["dat_selectstation"])

header = ["id", "res", "event", "stage 1", "stage 2", "stage 3"]
mdtable = Matrix{String}(undef, length(events), length(header))

# complex
# for ie = eachindex(events)
#     erdir = joinpath(PTH["dat_selectstation"], events[ie], "result")
#     mdtable[ie, 1] = string(ie)
#     mdtable[ie, 2] = string(length(events) - ie + 1)
#     mdtable[ie, 3] = events[ie]
#     for j = 1:3
#         fjld = joinpath(erdir, "result_stage$j.jld2")
#         ffig = relpath(joinpath(erdir, "result_stage$(j)_mt.png"))
#         if isfile(fjld)
#             if isfile(ffig)
#                 mdtable[ie, 3+j] = "[$(checksym)]($(ffig))"
#             else
#                 mdtable[ie, 3+j] = checksym
#             end
#         else
#             mdtable[ie, 3+j] = ""
#         end
#     end
# end

# simple
nrunning = 0
for ie = eachindex(events)
    global nrunning
    erdir = joinpath(PTH["dat_selectstation"], events[ie], "result")
    mdtable[ie, 1] = string(ie)
    mdtable[ie, 2] = string(length(events) - ie + 1)
    local flag = 0
    for j = 1:3
        fjld = joinpath(erdir, "result_stage$j.jld2")
        if isfile(fjld)
            mdtable[ie, 3+j] = checksym
            flag += 1
        else
            mdtable[ie, 3+j] = ""
        end
    end
    if events[ie] in errorevent
        mdtable[ie, 3] = "x " * events[ie]
        continue
    end
    if flag < 3 && isfile(joinpath(erdir, "begin.txt"))
        nrunning += 1
        mdtable[ie, 3] = "$(nrunning)% " * events[ie]
    else
        mdtable[ie, 3] = events[ie]
    end
end

function strlenwith(s::AbstractString, n::Integer)
    slen = length(s)
    if n <= slen
        return String(s[1:n])
    end
    m = floor(Int, (n - slen) / 2)
    return String(" "^m * s * " "^(n - m - slen))
end

function printmdline(io::IO, str::AbstractVector{<:AbstractString}, len::Vector{<:Integer})
    println(io, "|", join(map((_s, _i) -> strlenwith(_s, _i), str, len), '|'), "|")
end

setalign(n::Integer) = String(":" * ("-"^(n - 4)) * ":")

lenset = permutedims(header) |> _t->vcat(_t, mdtable) |> _t -> length.(_t) |> _t -> maximum(_t, dims=1) |> vec
lenset .+= 2

open("inverse_status.md", "w") do io
    printmdline(io, header, lenset)
    printmdline(io, setalign.(lenset), lenset)
    for l in eachrow(mdtable)
        printmdline(io, l, lenset)
    end
    println(io, "")
end
