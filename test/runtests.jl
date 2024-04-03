using DrWatson
@quickactivate "CSDNoise"

using TestItemRunner
using CSDNoise

println("Starting tests")
ti = time()

@run_package_tests filter = ti -> (occursin("CSDNoise", ti.filename)) verbose =
    true

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60; digits = 3), " minutes")
