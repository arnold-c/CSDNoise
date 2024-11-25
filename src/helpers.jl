using DrWatson
@quickactivate "CSDNoise"

outdir(args...) = projectdir("out", args...)

function sentencecase(s)
    @assert length(s) >= 2
    return uppercase(s[1]) * s[2:end]
end
