export
    real_sampler,
    complex_sampler

function real_sampler(n_samples::Int64,sample_length::Int64)
    [randn(Float64,sample_length) for i in 1:n_samples]
end

function complex_sampler(n_samples::Int64,sample_length::Int64)
    [randn(ComplexF64,sample_length) for i in 1:n_samples]
end