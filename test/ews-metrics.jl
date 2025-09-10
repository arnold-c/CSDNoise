@testset "ews-metrics.jl" begin
    using CSDNoise, StatsBase, Dates

    @testset "Timeseries aggregation" begin
        testvec = collect(1:10)
        @test isequal(aggregate_timeseries(testvec, Day(1)), testvec)
        @test isequal(
            aggregate_timeseries(testvec, Day(2)),
            [sum(1:2), sum(3:4), sum(5:6), sum(7:8), sum(9:10)],
        )
        @test isequal(
            aggregate_timeseries(testvec, Day(3)),
            [sum(1:3), sum(4:6), sum(7:9)],
        )
        @test isequal(
            aggregate_timeseries(testvec, Day(3), StatsBase.mean),
            [StatsBase.mean(1:3), StatsBase.mean(4:6), StatsBase.mean(7:9)],
        )

        @test isequal(
            aggregate_timeseries(testvec, Week(1)), [sum(1:7)]
        )
        @test isequal(
            aggregate_timeseries(testvec, Week(1), StatsBase.mean),
            [StatsBase.mean(1:7)],
        )
    end

    @testset "Mean functions" begin
        # Test Centered method branch coverage
        # Branch 1: i < bandwidth && i + bandwidth <= tlength (normal case)
        short_vec = [1.0, 2.0, 3.0, 4.0, 5.0]
        centered_short = spaero_mean(EWSMethod(Centered()), short_vec, 3)
        @test centered_short[1] ≈ 2.0  # mean([1,2,3]) - i=1, bandwidth=3, i+bandwidth=4 <= 5
        @test centered_short[2] ≈ 2.5  # mean([1,2,3,4]) - i=2, bandwidth=3, i+bandwidth=5 <= 5

        # Branch 2: i < bandwidth but i + bandwidth > tlength (short series)
        very_short_vec = [1.0, 2.0, 3.0, 4.0]
        centered_very_short = spaero_mean(EWSMethod(Centered()), very_short_vec, 3)
        @test centered_very_short[1] ≈ 2.0  # mean([1,2,3]) - i=1, bandwidth=3, i+bandwidth=4 <= 4
        @test centered_very_short[2] ≈ 2.5  # mean([1,2,3,4]) - i=2, bandwidth=3, i+bandwidth=5 > 4, uses entire series

        # Branch 3: i >= bandwidth but i + bandwidth > tlength (near end)
        @test centered_short[3] ≈ 3.0  # mean([1,2,3,4,5]) - i=3, bandwidth=3, i+bandwidth=6 > 5, uses from (3-3+1)=1 to end
        @test centered_short[4] ≈ 3.5  # mean([2,3,4,5]) - i=4, bandwidth=3, i+bandwidth=7 > 5, uses from (4-3+1)=2 to end
        @test centered_short[5] ≈ 4.0  # mean([3,4,5]) - i=5, bandwidth=3, i+bandwidth=8 > 5, uses from (5-3+1)=3 to end

        # Branch 4: else case (full centered window) - need longer series to trigger this
        longer_vec = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
        centered_long = spaero_mean(EWSMethod(Centered()), longer_vec, 3)
        @test centered_long[4] ≈ 4.0  # mean([2,3,4,5,6]) - i=4, bandwidth=3, full centered window from (4-3+1)=2 to (4+3-1)=6

        # Test that EWSMethod wrapper works correctly
        # The direct variant calls don't have methods - only EWSMethod does

        # Test Backward method branch coverage
        # Branch 1: i < bandwidth (partial window)
        backward_partial = spaero_mean(EWSMethod(Backward()), longer_vec, 4)
        @test backward_partial[1] ≈ 1.0  # mean([1])
        @test backward_partial[2] ≈ 1.5  # mean([1,2])
        @test backward_partial[3] ≈ 2.0  # mean([1,2,3])

        # Branch 2: i >= bandwidth (full backward window)
        @test backward_partial[4] ≈ 2.5  # mean([1,2,3,4])
        @test backward_partial[5] ≈ 3.5  # mean([2,3,4,5])
        @test backward_partial[6] ≈ 4.5  # mean([3,4,5,6])

        # Test edge case: bandwidth = 1
        mean_vec_bw1 = spaero_mean(EWSMethod(Backward()), longer_vec, 1)
        @test mean_vec_bw1 ≈ Float64.(longer_vec)  # Should equal original values

        # Test in-place version
        result_vec = zeros(Float64, length(longer_vec))
        spaero_mean!(result_vec, EWSMethod(Backward()), longer_vec, 2)
        mean_vec_backward = spaero_mean(EWSMethod(Backward()), longer_vec, 2)
        @test result_vec ≈ mean_vec_backward

        # Additional edge cases
        # Test with bandwidth equal to series length
        single_elem = [5.0]
        mean_single_backward = spaero_mean(EWSMethod(Backward()), single_elem, 1)
        @test mean_single_backward[1] ≈ 5.0

        mean_single_centered = spaero_mean(EWSMethod(Centered()), single_elem, 1)
        @test mean_single_centered[1] ≈ 5.0

        # Test with bandwidth larger than series length
        two_elem = [1.0, 3.0]
        mean_two_backward = spaero_mean(EWSMethod(Backward()), two_elem, 5)
        @test mean_two_backward[1] ≈ 1.0  # mean([1])
        @test mean_two_backward[2] ≈ 2.0  # mean([1,3])

        mean_two_centered = spaero_mean(EWSMethod(Centered()), two_elem, 5)
        @test mean_two_centered[1] ≈ 2.0  # mean([1,3]) - uses entire series
        @test mean_two_centered[2] ≈ 2.0  # mean([1,3]) - uses entire series

        # Test EWSMethod wrappers return vecs of correct length
        @test length(centered_long) == length(longer_vec)
        @test length(backward_partial) == length(longer_vec)
    end

    @testset "Backwards/Centered offset relationship" begin
        # For long enough time series, backwards and centered methods should have
        # identical values where their windows overlap exactly
        long_vec = collect(1.0:20.0)
        bandwidth = 5

        backward_result = spaero_mean(EWSMethod(Backward()), long_vec, bandwidth)
        centered_result = spaero_mean(EWSMethod(Centered()), long_vec, bandwidth)

        # The key insight: backward[bandwidth] uses window [1:bandwidth]
        # and centered[1] also uses window [1:bandwidth] (when 1 < bandwidth and 1+bandwidth <= length)
        # This is the only case where the windows are identical

        # Test the specific case where windows are identical:
        # backward[bandwidth] uses [1:bandwidth], centered[1] uses [1:bandwidth]
        @test backward_result[bandwidth] ≈ centered_result[1] atol = 1.0e-10

        # Test with a simple example to verify our understanding
        simple_vec = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        simple_bandwidth = 3

        simple_backward = spaero_mean(EWSMethod(Backward()), simple_vec, simple_bandwidth)
        simple_centered = spaero_mean(EWSMethod(Centered()), simple_vec, simple_bandwidth)

        # backward[3] uses [1,2,3], centered[1] uses [1,2,3] -> should be equal
        @test simple_backward[simple_bandwidth] ≈ simple_centered[1] atol = 1.0e-10

        # More generally, test that there exists a systematic offset relationship
        # in the stable regions of both methods

        # For longer series, test that the backwards method and the centered method are related
        # when the backwards method add 1 to the bandwidth (to account for the central value in
        # the centered method) and the centered method halves the bandwidth (to account for the
        # window being constructed with the bandwidth on either side of the central value)
        very_long_vec = collect(1.0:50.0)
        large_bandwidth = 10

        very_long_backward = spaero_mean(EWSMethod(Backward()), very_long_vec, large_bandwidth + 1)
        very_long_centered = spaero_mean(EWSMethod(Centered()), very_long_vec, large_bandwidth ÷ 2)

        stable_backward = very_long_backward[(large_bandwidth + 1):(end)]
        stable_centered = very_long_centered[(large_bandwidth ÷ 2 + 1):(end - large_bandwidth ÷ 2)]

        @test stable_backward == stable_centered
    end


    @testset "Variance functions" begin
        testvec = collect(1:10)

        # Test variance calculation with pre-computed mean
        mean_vec = spaero_mean(EWSMethod(Backward()), testvec, 3)
        var_vec = spaero_var(EWSMethod(Backward()), mean_vec, testvec, 3)
        @test length(var_vec) == length(testvec)
        @test all(var_vec .>= 0)  # Variance should be non-negative

        # Test variance calculation without pre-computed mean
        var_vec2 = spaero_var(EWSMethod(Backward()), testvec, 3)
        @test var_vec ≈ var_vec2

        # Test with constant vector (variance should be 0)
        const_vec = fill(5.0, 10)
        var_const = spaero_var(EWSMethod(Backward()), const_vec, 3)
        @test all(var_const .≈ 0.0)
    end

    @testset "Higher moments" begin
        testvec = collect(1:10)

        # Test skewness
        skew_vec = spaero_skew(EWSMethod(Backward()), testvec, 3)
        @test length(skew_vec) == length(testvec)

        # Test kurtosis
        kurt_vec = spaero_kurtosis(EWSMethod(Backward()), testvec, 3)
        @test length(kurt_vec) == length(testvec)

        # Test helper functions with known values
        m3_vec = fill(8.0, 5)  # Third moment
        sd3_vec = fill(2.0, 5)  # Standard deviation cubed
        skew_result = spaero_skew(m3_vec, sd3_vec)
        @test all(skew_result .≈ 4.0)  # 8/2 = 4

        m4_vec = fill(16.0, 5)  # Fourth moment
        var2_vec = fill(4.0, 5)  # Variance squared
        kurt_result = spaero_kurtosis(m4_vec, var2_vec)
        @test all(kurt_result .≈ 4.0)  # 16/4 = 4
    end

    @testset "CoV and IoD" begin
        # Test coefficient of variation
        sd_vec = [1.0, 2.0, 3.0]
        mean_vec = [2.0, 4.0, 6.0]
        cov_result = spaero_cov(sd_vec, mean_vec)
        @test cov_result ≈ [0.5, 0.5, 0.5]  # sd/mean

        # Test index of dispersion
        var_vec = [1.0, 4.0, 9.0]
        iod_result = spaero_iod(var_vec, mean_vec)
        @test iod_result ≈ [0.5, 1.0, 1.5]  # var/mean

        # Test full calculation
        testvec = [1.0, 2.0, 3.0, 4.0, 5.0]
        cov_full = spaero_cov(EWSMethod(Backward()), testvec, 2)
        @test length(cov_full) == length(testvec)

        iod_full = spaero_iod(EWSMethod(Backward()), testvec, 2)
        @test length(iod_full) == length(testvec)
    end

    @testset "Autocorrelation and Autocovariance" begin
        testvec = collect(1:10)

        # Test autocovariance with lag 1
        autocov_vec = spaero_autocov(EWSMethod(Backward()), testvec, 3; lag = 1)
        @test length(autocov_vec) == length(testvec)
        @test isnan(autocov_vec[1])  # First lag should be NaN

        # Test with pre-computed mean
        mean_vec = spaero_mean(EWSMethod(Backward()), testvec, 3)
        autocov_vec2 = spaero_autocov(EWSMethod(Backward()), mean_vec, testvec, 3; lag = 1)
        # Compare non-NaN values since NaN ≈ NaN is false
        @test autocov_vec[2:end] ≈ autocov_vec2[2:end]
        @test isnan(autocov_vec[1]) && isnan(autocov_vec2[1])

        # Test autocorrelation
        autocor_vec = spaero_autocor(EWSMethod(Backward()), testvec, 3; lag = 1)
        @test length(autocor_vec) == length(testvec)
        @test isnan(autocor_vec[1])  # First lag should be NaN

        # Test autocorrelation with pre-computed values
        sd_vec = sqrt.(spaero_var(EWSMethod(Backward()), testvec, 3))
        autocor_vec2 = spaero_autocor(autocov_vec, sd_vec; lag = 1)
        # Test precomputed and wrapper functions produce the same values
        # First two are NaN so ignore
        @test autocor_vec[3:end] ≈ autocor_vec2[3:end]

        # Test with different lags
        autocov_lag2 = spaero_autocov(EWSMethod(Backward()), testvec, 3; lag = 2)
        @test isnan(autocov_lag2[1]) && isnan(autocov_lag2[2])  # First two should be NaN
    end

    @testset "Kendall correlation" begin
        # Test with increasing trend
        increasing_vec = [1.0, 2.0, 3.0, 4.0, 5.0]
        tau_inc = spaero_corkendall(increasing_vec)
        @test tau_inc > 0  # Should be positive correlation

        # Test with decreasing trend
        decreasing_vec = [5.0, 4.0, 3.0, 2.0, 1.0]
        tau_dec = spaero_corkendall(decreasing_vec)
        @test tau_dec < 0  # Should be negative correlation

        # Test with NaN values (should be filtered)
        vec_with_nan = [1.0, NaN, 3.0, 4.0, 5.0]
        tau_nan = spaero_corkendall(vec_with_nan)
        @test !isnan(tau_nan)  # Should not be NaN after filtering
        @test tau_nan > 0

        # Test with constant values
        constant_vec = fill(3.0, 5)
        tau_const = spaero_corkendall(constant_vec)
        @test isnan(tau_const)
    end

    @testset "Bandwidth aggregation" begin
        using CSDNoise: EWSMetricSpecification

        # Test basic aggregation
        spec = EWSMetricSpecification(
            EWSMethod(Backward()),
            Day(2),  # aggregation
            Day(10), # bandwidth
            1  # lag
        )
        agg_bw = CSDNoise.aggregate_bandwidth(spec)
        @test agg_bw == 5  # 10 ÷ 2

        # Test with Day(1) (no aggregation)
        spec2 = EWSMetricSpecification(
            EWSMethod(Backward()),
            Day(1),  # aggregation
            Day(6),  # bandwidth
            1
        )
        agg_bw2 = CSDNoise.aggregate_bandwidth(spec2)
        @test agg_bw2 == 6

        # Test assertion for invalid bandwidth
        spec_invalid = EWSMetricSpecification(
            EWSMethod(Backward()),
            Day(4),  # aggregation > bandwidth
            Day(3),  # bandwidth
            1
        )
        @test_throws AssertionError CSDNoise.aggregate_bandwidth(spec_invalid)
    end
end
