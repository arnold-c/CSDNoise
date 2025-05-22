default: tests manuscript

manuscript:
	julia ./manuscript/scripts/optimal-thresholds.jl
	typst compile ./manuscript/combined-manuscript.typ

tests:
	julia --project=. -e "using Pkg; Pkg.test()"
