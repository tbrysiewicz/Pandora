#export Sampler

abstract type Sampler end

function real_uniform(n::Int,n_samples::Int)
    [rand(Float64,n) for i in 1:n_samples]
end

function real_normal(n::Int,n_samples::Int)
    [randn(Float64,n) for i in 1:n_samples]
end

function complex_uniform(n::Int,n_samples::Int)
[rand(ComplexF64,n) for i in 1:n_samples]
end

function complex_normal(n::Int,n_samples::Int)
[randn(ComplexF64,n) for i in 1:n_samples]
end

mutable struct EllipseSampler{T<: Number} <: Sampler
    n :: Int                                            # n is the dimension of the parameter space, i.e., the length of a vector in the sample produced. 
    n_samples :: Int
    predistribution 
    transform_matrix :: Matrix{T} 
    translation :: Vector{T}
end

function n(ell_sampler::EllipseSampler)
    ell_sampler.n
end

function n_samples(ell_sampler::EllipseSampler)
    ell_sampler.n_samples
end

function predistribution(ell_sampler::EllipseSampler)
    ell_sampler.predistribution
end

function transform_matrix(ell_sampler::EllipseSampler)
    ell_sampler.transform_matrix
end

function translation(ell_sampler::EllipseSampler)
    ell_sampler.translation
end

function sample(S::EllipseSampler)
    map(x->Vector{ComplexF64}(transform_matrix(S)*x+translation(S)), predistribution(S)(S.n,n_samples(S)))
end

function is_real(S::Sampler)
    s = sample(S)
    return(is_real(s[1]))
end


function initialize_real_sampler(n::Int,n_samples::Int;uniform=false)
    if uniform
        EllipseSampler(n,n_samples, real_uniform,Matrix{Float64}(I,n,n), zeros(n))
    else
        EllipseSampler(n,n_samples, real_normal,Matrix{Float64}(I,n,n), zeros(n))
    end
end


function initialize_complex_sampler(n::Int,n_samples::Int;uniform=false)
    if uniform
        EllipseSampler(n,n_samples, complex_uniform,Matrix{ComplexF64}(I,n,n), Vector{ComplexF64}(zeros(n)))
    else
        EllipseSampler(n,n_samples, complex_normal,Matrix{ComplexF64}(I,n,n), Vector{ComplexF64}(zeros(n)))
    end
end