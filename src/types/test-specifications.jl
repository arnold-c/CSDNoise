export IndividualTestSpecification

Base.@kwdef struct IndividualTestSpecification
    sensitivity::Float64
    specificity::Float64
    test_result_lag::Int64
end
