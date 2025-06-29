export 
    Sampler, 
    UnitSampler, 
    UniformSampler, 
    NormalSampler,
    field,
    AffineTransformation,
    TransformedSampler,
    translate,
    translate!,
    set_translation,
    set_transform_matrix,
    scale,
    scale!



abstract type Sampler end


function Base.show(io::IO, S::Sampler) 
    print(io, name(S))
end

"""
    UnitSampler{T}(n)

A sampler that produces an `n`-tuples of elements of type `T` via `rand(T,n)`
"""
struct UnitSampler{T<: Number} <: Sampler
    dimen :: Int                                        # n is the dimension of the parameter space, 
                                                        # i.e., the length of a vector in the sample produced.                                
end

name(US::UnitSampler{T}) where {T<:Number} = string("Unit sampler of unit box centered at origin in ", T, "^", US.dimen)


"""
    UniformSampler{T}(n)

    A sampler which produces a uniform distrubtion on the unit ball in `T^n`. 
"""
struct UniformSampler{T<: Number} <: Sampler
    dimen :: Int                                        # n is the dimension of the parameter space, 
                                                        # i.e., the length of a vector in the sample produced.                                
end
name(UNS::UniformSampler{T}) where {T<:Number} = string("Uniform sampler of unit ball in ", T, "^", UNS.dimen)


"""
    NormalSampler{T}(n)

A sampler that produces an `n`-tuples of elements of type `T` via `randn(T,n)`
"""
struct NormalSampler{T<: Number} <: Sampler
    dimen :: Int                                        # n is the dimension of the parameter space, 
                                                        # i.e., the length of a vector in the sample produced.                                
end

name(NS::NormalSampler{T}) where {T<:Number} = string("Normal sampler of unit box centered at origin in ", T, "^", NS.dimen)


field(::UnitSampler{T}) where T<: Number = T
field(::UniformSampler{T}) where T<: Number = T
field(::NormalSampler{T}) where T<: Number = T


function (UNS::UniformSampler{T})(n::Int) where {T<:Number}
    X = [randn(T,UNS.dimen) for _ in 1:n]
    X = [x./norm(x) for x in X]
    return [rand()^(1/UNS.dimen).*x for x in X]
end

function (US::UnitSampler{T})(n::Int) where {T<:Number}
    [[(-1)^rand([0,1])*r for r in rand(T,US.dimen)] for i in 1:n]
end

function (NS::NormalSampler{T})(n::Int) where {T<:Number}
    [randn(T,NS.dimen) for i in 1:n]
end

"""
    AffineTransformation(M::AbstractMatrix, t::AbstractVector)
Creates an affine transformation defined by a matrix `M` and a translation vector `t`.
"""
mutable struct AffineTransformation{T<: Number}
    transform_matrix :: Matrix{T}                       # The matrix that linear transforms our predistribution
                                                        # to an ellipse sample.
    translation :: Vector{T}                            # Translation applied to predistribution to get the 
                                                        # desired ellipse sample.
end

function AffineTransformation(M::AbstractMatrix, t::AbstractVector)
    T = promote_type(eltype(M), eltype(t))
    newM = convert(Matrix{T}, M)
    newt = convert(Vector{T}, t)
    return AffineTransformation{T}(newM, newt)
end

"""
    TransformedSampler{T}(predistribution :: Sampler , affine_transformation :: AffineTransformation{T})

A sampler that applies an affine transformation to samples generated by a predistribution sampler.
"""
mutable struct TransformedSampler{T<: Number} <: Sampler
    predistribution :: Sampler                          # predistribution is the method used to sample points
                                                        # from a circle, before it is changed to an ellipse in
                                                        # our desired directions
    affine_transformation :: AffineTransformation{T}   # The affine transformation applied to the predistribution
end

name(TS::TransformedSampler{T}) where {T<:Number} = string("Transformation of (", name(predistribution(TS)), ") by linear transformation \n", TS.affine_transformation.transform_matrix, "\n and translation \n", TS.affine_transformation.translation)

predistribution(TS::TransformedSampler{T}) where {T<:Number} = TS.predistribution

function (TS::TransformedSampler{T})(n::Int) where {T<:Number}
    # Sample from the predistribution
    samples = predistribution(TS)(n)
    
    # Apply the affine transformation to each sample
    return [TS.affine_transformation.transform_matrix * x + TS.affine_transformation.translation for x in samples]
end

field(TS::TransformedSampler{T}) where {T<:Number} = T


function scale(TS::TransformedSampler{T}, s::Number) where {T<:Number}
    newM = TS.affine_transformation.transform_matrix * s
    return TransformedSampler{T}(predistribution(TS), AffineTransformation(newM, TS.affine_transformation.translation))
end
function scale!(TS::TransformedSampler{T}, s::Number) where {T<:Number}
    TS.affine_transformation.transform_matrix = TS.affine_transformation.transform_matrix * s
end

function set_translation(TS::TransformedSampler{T}, t::AbstractVector) where {T<:Number}
    newt = convert(Vector{T}, t)
    return TransformedSampler{T}(predistribution(TS), AffineTransformation(TS.affine_transformation.transform_matrix, newt))
end

function set_translation!(TS::TransformedSampler{T}, t::AbstractVector) where {T<:Number}
    newt = convert(Vector{T}, t)
    TS.affine_transformation.translation = newt
end

function set_transform_matrix(TS::TransformedSampler{T}, M::AbstractMatrix) where {T<:Number}
    newM = convert(Matrix{T}, M)
    return TransformedSampler{T}(predistribution(TS), AffineTransformation(newM, TS.affine_transformation.translation))
end

function set_transform_matrix!(TS::TransformedSampler{T}, M::AbstractMatrix) where {T<:Number}
    newM = convert(Matrix{T}, M)
    TS.affine_transformation.transform_matrix = newM
end


function translate(S::Sampler, t::AbstractVector)
    T = field(S)
    newt = convert(Vector{T}, t)
    if S isa TransformedSampler
        return TransformedSampler{T}(predistribution(S), AffineTransformation(S.affine_transformation.transform_matrix, S.affine_transformation.translation + newt))
    elseif S isa UnitSampler || S isa UniformSampler || S isa NormalSampler
        return TransformedSampler{T}(S, AffineTransformation(Matrix{Float64}(I, S.dimen, S.dimen), newt))
    else
        throw(ArgumentError("Cannot translate sampler of type $(typeof(S))"))
    end
end

function translate!(TS::TransformedSampler{T}, t::AbstractVector) where {T<:Number}
    newt = convert(Vector{T}, t)
    TS.affine_transformation.translation += newt
end

Base.convert(::Type{TransformedSampler}, S::Sampler)  = translate(S, zeros(field(S), S.dimen)) 