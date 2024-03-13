using DrWatson, TestItemRunner
@quickactivate "CSDNoise"

using CSDNoise

# Here you include files using `srcdir`
# include(srcdir("file.jl"))

# Run test suite
println("Starting tests")
ti = time()

@run_package_tests verbose = true

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60; digits = 3), " minutes")
