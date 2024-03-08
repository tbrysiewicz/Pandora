#TODO: Get progress meter working 
#TODO: Make sure proper attention to irreducibility is given
#TODO: Make it so that one can construct the enumerative problems associated to bases

export
    algebraic_matroid,
    affine_multidegree,
    affine_multidegrees,
    condition_number_of_basis,
    condition_numbers_of_candidate_bases,
    numerical_algebraic_matroid,
    numerical_bases

function condition_number_of_basis(V :: Variety, c :: Vector{Int}, p :: Vector{ComplexF64}, n :: Int)
    J = HomotopyContinuation.jacobian(system(V), p)
    C = LinearAlgebra.cond(J[:,filter(x->in(x,c)==false,1:n)])
    return(C)
end

function condition_numbers_of_candidate_bases(V :: Variety; dim = nothing)
    if is_populated(V)==false 
        if dim != nothing
            populate_one_point!(V,dim)
        else
            populate_witness!(V)
        end
    end
    n = ambient_dimension(V)
    dimension = Pandora.dim(V)
    candidateBases = Combinatorics.combinations(1:n, dimension)
    p = witness_points(V)[1]
    D = Dict{Vector{Int},Float64}()
    counter=0

    J = HomotopyContinuation.jacobian(system(V), p)

    for c in candidateBases
        counter=counter+1
        if floor(counter/1000)==counter//1000
            println(counter)
        end
        D[c]=LinearAlgebra.cond(J[:,filter(x->in(x,c)==false,1:n)])
    end
    return(D)
end

function numerical_bases(V :: Variety; tol = 1e7)
    D = condition_numbers_of_candidate_bases(V)
    bases = []
    for k in keys(D)
        if D[k]<tol
            push!(bases,k)
        end
    end
    return(Vector{Vector{Int}}(bases))
end

@doc raw"""
    numerical_algebraic_matroid(V::Variety; tol = 1e7)

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
function numerical_algebraic_matroid(V :: Variety; tol = 1e7)
    Bases = numerical_bases(V;tol=tol)
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
function affine_multidegree1(W::WitnessSet,I::Vector{Int64})
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
