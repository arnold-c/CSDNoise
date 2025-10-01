export outdir,
    sentencecase

outdir(args...) = DrWatson.projectdir("out", args...)

function sentencecase(s)
    @assert length(s) >= 2
    return uppercase(s[1]) * s[2:end]
end
