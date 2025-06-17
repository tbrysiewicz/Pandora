export Sampler, EllipseSampler

abstract type Sampler end

# Getter for Sampler

function n_samples(S::Sampler) 
    S.n_samples
end

# A concrete subtype of the type Sampler, that samples points 
# from an ellipse, say for example, by using an ellipse whose 
# axes are directed matching to the the previous step in sampling.

mutable struct EllipseSampler{T<: Number} <: Sampler
    n :: Int                                            # n is the dimension of the parameter space, 
                                                        # i.e., the length of a vector in the sample produced. 
    n_samples :: Int                                    # n_samples is the number of sample points.
    predistribution                                     # predistribution is the method used to sample points
                                                        # from a circle, before it is changed to an ellipse in
                                                        # our desired directions
    transform_matrix :: Matrix{T}                       # The matrix that linear transforms our predistribution
                                                        # to an ellipse sample.
    translation :: Vector{T}                            # Translation applied to predistribution to get the 
                                                        # desired ellipse sample.
end

# Base.show for ScoringScheme

function Base.show(io::IO, S::EllipseSampler)
    tenspaces="          "
    print(io,"\n")
    println(io, tenspaces, "A sampler with the following values  ")
    println("---------------------------------------------------------")
    println(io,"Dimesion of the parameter space    : ", S.n )
    println(io,"Number of samples                  : ", S.n_samples)
    println(io, "Predistribution                    : ", S.predistribution)
    println(io, "\nTransform_matrix: \n")
    display(S.transform_matrix)
    println(io, "\nTranslation vector: \n")
    display(S.translation)
end


#Getters for EllipseSampler

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

# Function that transforms our predistribution to our required
# ellipse sample by using the fields in S::EllipseSampler.

function sample(S::EllipseSampler)
    map(x->Vector{ComplexF64}(transform_matrix(S)*x+translation(S)), predistribution(S)(S.n,n_samples(S)))
end


# examples for predistributions in an EllipseSampler.

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

# checks whether the ellipse sample obtained is real.

function is_real(S::Sampler)
    s = sample(S)
    return(isreal(s[1]))
end

# Constructors for the type EllipseSampler:
import LinearAlgebra
I = LinearAlgebra.I
# For the real case i.e., when we want the predistribution to be real.
function initialize_real_sampler(n::Int,n_samples::Int;uniform=false)
    if uniform
        EllipseSampler(n,n_samples, real_uniform,Matrix{Float64}(I,n,n), zeros(n))
    else
        EllipseSampler(n,n_samples, real_normal,Matrix{Float64}(I,n,n), zeros(n))
    end
end

# For the complex case, i.e., when we want the predistribution to be complex.
function initialize_complex_sampler(n::Int,n_samples::Int;uniform=false)
    if uniform
        EllipseSampler(n,n_samples, complex_uniform,Matrix{ComplexF64}(I,n,n), Vector{ComplexF64}(zeros(n)))
    else
        EllipseSampler(n,n_samples, complex_normal,Matrix{ComplexF64}(I,n,n), Vector{ComplexF64}(zeros(n)))
    end
end

