#TODO: Get progress meter working 
#TODO: Make sure proper attention to irreducibility is given
#TODO: Make it so that one can construct the enumerative problems associated to bases

include("HelperFunctions.jl")


export
    algebraic_matroid,
    affine_multidegree,
    affine_multidegrees,
    condition_number_of_basis,
    condition_numbers_of_candidate_bases,
    numerical_algebraic_matroid,
    numerical_bases



mutable struct condBasis
    basis :: Vector{Int}
    conditionNum :: Float64
end

function Base.show(io::IO, cb::condBasis)
    print(cb.basis," -> ",cb.conditionNum)
end


function condition_numbers_of_candidate_bases(V :: Variety, jacobian :: Matrix{ComplexF64}, dim = nothing)

    ambientDim :: Int = ambient_dimension(V)

    conditionNums = Array{condBasis}(undef, binomial(ambientDim, dim))

    groundSet = collect(1:ambientDim)
    numCandidates = binomial(ambientDim, dim)

    Threads.@threads for i in ProgressBar(1:numCandidates)

        c :: Vector{Int} = rank_r_combination(ambientDim, dim, i)
            
        num :: Float64 = LinearAlgebra.cond(jacobian[:,setdiff(groundSet, c)])

        conditionNums[i] = condBasis(c,num)

    end

    return(conditionNums)

end


function numerical_bases(V :: Variety; dimension = nothing, amplify = 1)

    # first get the points that we need depending on if the dimension is given or not
    if is_populated(V) && length(witness_points(V))>=amplify# best case scenario the points already exist
        witnessPoints = witness_points(V)
        dimension = Pandora.dim(V)
    elseif is_populated(V) || (!isnothing(dimension)) #If we are populated, but not enough witness points (or dimension is given)
        witnessPoints = get_n_points(V, dimension, amplify)
    else # worst case we don't know points or dimension so we have to find the whole witness set
        witnessPoints = witness_points(V)
        dimension = Pandora.dim(V)
    end

    # make the jacobian and pass it off for the first run of condition numbers
    jac :: Matrix{ComplexF64} = HomotopyContinuation.jacobian(system(V), witnessPoints[1])
    conditionNums = condition_numbers_of_candidate_bases(V, jac, dimension)

    # put the condition numbers in a matrix so we can cluster them and find the tolerence
    conditionNumsMatrix :: Matrix{Float64} = reshape([((x) -> isfinite(x) ? log10(x) : 308.0)(c.conditionNum) for c in conditionNums], 1, length(conditionNums))
    clusters = kmeans(conditionNumsMatrix, 2)
    tolerence :: Float64 = 10 ^ ((clusters.centers[1] + clusters.centers[2])/2)

    # now we can execute the amplification procedure
    for i in 2:amplify
        # first create the new jacobian with the new point
        jac = HomotopyContinuation.jacobian(system(V), witnessPoints[i])
        
        ambientDim :: Int = ambient_dimension(V)    
        groundSet = collect(1:ambientDim)

        # look at the condition nums
        for i in conditionNums
            if i.conditionNum > tolerence
                newCond = LinearAlgebra.cond(jac[:,setdiff(groundSet, i.basis)])
                #Consider this?
                #m= min(newCond,i.conditionNum)
                #i.conditionNum = m
                if (i.conditionNum - newCond) > 0.0
                    i.conditionNum = newCond
                end

            end
        end
        
        # once large condition nums have been rechecked redo the clustering
        conditionNumsMatrix = reshape([((x) -> isfinite(x) ? log10(x) : 308.0)(c.conditionNum) for c in conditionNums], 1, length(conditionNums))
        clusters = kmeans(conditionNumsMatrix, 2)
        tolerence = 10 ^ ((clusters.centers[1] + clusters.centers[2])/2)

    end

    bases :: Vector{Vector{Int}} = []

    for k in conditionNums

        if k.conditionNum < tolerence
            push!(bases,k.basis)
        end

    end

    return(bases)

end


@doc raw"""
    numerical_algebraic_matroid(V::Variety)

 Returns the algebraic matroid of the variety V computed numerically
 # Examples
 ```jldoctest
 julia> V = SpecialOrthogonalGroup(3);

 julia> M = numerical_algebraic_matroid(V)
 Matroid of rank 3 on 9 elements
 
 julia> nonbases(M)
 6-element Vector{Vector{Int64}}:
  [1, 2, 3]
  [1, 4, 7]
  [2, 5, 8]
  [3, 6, 9]
  [4, 5, 6]
  [7, 8, 9]
 ```
 """
function numerical_algebraic_matroid(V :: Variety)

    Bases = numerical_bases(V)

    M = Oscar.matroid_from_bases(Bases,ambient_dimension(V))

    return(M)

end

#TODO: make this a symbolic computation (find exact symbolic det function and use GB for zero-test (ideal containment))
@doc raw"""
    algebraic_matroid(V::Variety)

 Returns the algebraic matroid of the variety V
 # Examples
 ```jldoctest
 julia> V = SpecialOrthogonalGroup(3);

 julia> M = algebraic_matroid(V)
 Matroid of rank 3 on 9 elements
 
 julia> nonbases(M)
 6-element Vector{Vector{Int64}}:
  [1, 2, 3]
  [1, 4, 7]
  [2, 5, 8]
  [3, 6, 9]
  [4, 5, 6]
  [7, 8, 9]
 ```
 """
function algebraic_matroid(V::Variety)
    AM = affine_multidegrees(V)
    Vars = variables(system(V))
    Bases = [findall(x->x ∈ c,Vars) for c in keys(AM)]
    Oscar.matroid_from_bases(Bases,ambient_dimension(V))
end


#TODO: Update these so that the bases are calculated fast, and monodromy is used to find degrees?
# Possibly monodromy recover would be an optional argument and one can do the big computation if
# reliability is suspect for monodromy recover.
@doc raw"""
    affine_multidegrees(V::Variety)

 Returns a dictionary from the bases of the algebraic matroid of $V$ to the degree of the corresponding branched cover.
 # Examples
 ```jldoctest
 julia> V = SpecialOrthogonalGroup(3);

 julia> affine_multidegrees(V)
Dict{Vector{Variable}, Int64} with 78 entries:
  [x₂₋₁, x₁₋₃, x₃₋₃] => 4
  [x₃₋₁, x₂₋₂, x₂₋₃] => 4
  [x₁₋₁, x₁₋₂, x₃₋₂] => 4
  [x₂₋₁, x₁₋₂, x₁₋₃] => 4
  [x₁₋₁, x₃₋₁, x₃₋₂] => 4
  [x₁₋₁, x₂₋₂, x₁₋₃] => 4
  [x₁₋₁, x₃₋₁, x₁₋₂] => 4
  [x₁₋₂, x₃₋₂, x₂₋₃] => 4
  [x₃₋₁, x₁₋₃, x₂₋₃] => 4
  [x₃₋₁, x₂₋₃, x₃₋₃] => 4
  [x₁₋₁, x₁₋₃, x₃₋₃] => 4
  [x₁₋₁, x₂₋₁, x₂₋₃] => 4
  [x₁₋₁, x₂₋₁, x₃₋₃] => 4
  [x₁₋₁, x₃₋₂, x₂₋₃] => 8
  [x₁₋₁, x₂₋₂, x₃₋₃] => 8
  ⋮                  => ⋮
 ```
 """
function affine_multidegrees(V::Variety)
    F = system(V)
    if is_populated(V)==false
        populate_witness!(V)
    end
    Vars = variables(F)
    n = length(Vars)
    d = dim(V)
    C = collect(combinations(1:n,d))
    AM = Dict{Vector{Variable},Int64}()
    W = witness_set(V)
#    println("Checking ",length(C)," many potential multidegrees.")
    counter=0
    talker=0
    printer_count=length(C)/1000
    for c in C
        counter=counter+1
        talker=talker+1
        if talker>printer_count
  #          println(counter,"/",length(C))
            talker=0
        end
        amd=affine_multidegree(W,c)
        if amd!=0
            AM[Vars[c]]=amd
        end
    end
    return(AM)
end




@doc raw"""
    affine_multidegree(W::WitnessSet,I::Vector{Int64}))

 Returns the degree of the branched cover from V to the variables indexed by I.
 # Examples
 ```jldoctest
 julia> V = SpecialOrthogonalGroup(3);

 julia> affine_multidegree(witness_set(V),[1,2,4])
 4
 ```
 """
function affine_multidegree(W::WitnessSet,I::Vector{Int64})
    d = HomotopyContinuation.dim(W) #Dimension of witness set
    n = length(variables(W.F)) #Ambient Dimension
    A = randn(ComplexF64,d,n) #The size of linear space on needs to slice with
    b = randn(ComplexF64,d) #The affine part of the linear space
    pt = [W.R[1]] 
    W2 = WitnessSet(W.F,W.L,pt)
    @assert length(I) == d
    for i in 1:d
        for j in 1:n
            if I[i]==j
                A[i,j]=1.0
            else
                A[i,j]=0.0
            end
        end
    end
    b = A*solution(pt[1]) + 0.0001*randn(ComplexF64,d)
    L = LinearSubspace(A,b)
    SpecialW = witness_set(W,LinearSubspace(A,b)) #solving
    MS = monodromy_solve(SpecialW.F, solutions(SpecialW.R), SpecialW.L)
    return(length(solutions(MS)))

end

