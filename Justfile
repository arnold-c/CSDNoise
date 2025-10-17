default: optimizations

optimizations:
	julia --threads=auto --startup-file=no --history-file=no --project=CSDNoiseCore ./scripts/gridsearch-optimization.jl

manuscript:
	julia ./manuscript/scripts/optimal-thresholds.jl
	typst compile ./manuscript/combined-manuscript.typ

tests:
	julia --project=. -e "using Pkg; Pkg.test()"
