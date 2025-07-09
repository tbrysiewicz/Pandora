import Base: length, show, getindex
export 
    MonodromyLoop,
    Loop,
    large_monodromy_sample,
    monodromy_group,
    galois_group,
    monodromy_sample,
    permutation,
    permutation!,
    group_generated_by_monodromy_loops,
    is_decomposable, 
    orbits

"""
     Loop(P::Vector{Vector{ComplexF64}})

Creates a loop object from a vector of parameter vectors `P`. The first and last elements of `P` must be the same, representing the start and end of the loop.
"""
struct Loop
    P :: Vector{Vector{ComplexF64}} # A list of parameters starting and ending at F[2]
    function Loop(P::Vector{Vector{ComplexF64}})
        @assert length(P) >= 2 "Loop must have at least two parameters."
        @assert isapprox(P[1], P[end]) "Loop must start and end at the same parameter."
        new(P)
    end
end

function length(L::Loop)
    return length(L.P)
end

function getindex(L::Loop, i::Int)
    return L.P[i]
end


"""
    MonodromyLoop(F::Fibre, L::Loop, sigma::Union{PermGroupElem,Nothing})
"""
mutable struct MonodromyLoop
    F :: Fibre                      # Gives sols and ordering on them. (S,P)
    L :: Loop          # A list of parameters starting and ending at F[2]
    sigma :: Union{PermGroupElem,Nothing}

    function MonodromyLoop(F::Fibre, P::Vector{Vector{ComplexF64}}, sigma::Union{PermGroupElem,Nothing})
        @assert length(P) >= 2 "Monodromy loop must have at least two parameters."
        @assert isapprox(P[1], P[end]) "Monodromy loop must start and end at the same parameter."
        new(F, Loop(P), sigma)
    end

    function MonodromyLoop(F::Fibre, P::Vector{Vector{ComplexF64}})
        return MonodromyLoop(F, P, nothing)
    end

    function MonodromyLoop(EP::EnumerativeProblem, P::Vector{Vector{ComplexF64}})
        F = base_fibre(EP)
        
        if !isapprox(P[1],F[2])
            pushfirst!(P, F[2]) # Inconsistent loop: starting at given fibre parameter
        end
        if !isapprox(P[end],F[2])
            push!(P, F[2]) # Incomplete loop: using start parameter as end parameter
        end
        return MonodromyLoop(F, P, nothing)
    end
end

"""
    base_fibre(ML::MonodromyLoop)

    Returns the base fibre of the monodromy loop `ML`.
"""
function base_fibre(ML::MonodromyLoop)
    return ML.F
end

function Base.show(io::IO, ML::MonodromyLoop)
    if ML.sigma === nothing
        print(io,"A monodromy ",length(ML.L),"-loop.")
    else
        print(io,"A monodromy ",length(ML.L),"-loop with permutation ",ML.sigma)
    end
end


"""
    is_permutation(p::Vector{Int},d::Int64) :: Bool

    This function checks if the vector `p` represents a valid permutation of the integers from 1 to `d`.
    A valid permutation is injective and surjective, meaning it contains each integer from 1 to `d` exactly once.
    If `p` contains `nothing`, it indicates that the mapping is not surjective.
"""
function is_permutation(p::Union{Vector{Int},Vector{Union{Nothing,Int64}}},d::Int64) :: Bool
    u=unique(p) 
    return(!(length(u)<d || in(nothing,u))) #nothing indicates not surjective, length indicates not injective
end

"""
    numerical_function(S1::Vector{T}, S2::Vector{T}; tol::Float64 = 1e-6)

    This function determines the function `f: S2 -> S1` such that `f(s) = i` if `isapprox(S1[i], s, atol=tol)` and 
        represents it as a vector of indices. For example, the function `f: {a,b,c} -> {a,b,c}` for which `f(a)=a,f(b)=a,f(c)=b` would be represented as [1,1,2].
    If `f` represents a permutation, the output of this function is the one-line notation of the permutation.

    It is assumed that `S1` consists of distinct elements, and that `S2` is a (numerical) subset of `S1`.
    If `S2` is not a subset of `S1`, then certain entries of this vector will be `nothing`, indicating that the corresponding element of `S2` is not found in `S1`.
"""
function numerical_function(S1::Vector{T} where T, S2::Vector{T} where T; tol::Float64 = 1e-6) 
    one_line = [findfirst(x->isapprox(x,s),S1) for s in S2]
    return(one_line)
end

"""
    (EP::EnumerativeProblem)(ML::MonodromyLoop)

    Tracks the base fibre of the monodromy loop `ML` over the loop defined by `ML.L` with respect to the enumerative problem `EP`.
"""
function (EP::EnumerativeProblem)(ML::MonodromyLoop)
    (S,P) = base_fibre(ML)
    loop = ML.L.P
    for l in loop[2:end]
        (S,P) = (EP((S,P),l),l)
    end
    return(S)
end

"""
    (EP::EnumerativeProblem)(L::Loop)

    Tracks the base fibre `F = (S,P)` of `EP` over the loop  `L`, conjugated with a path from `P` to `L[1]=L[end]`.
"""
function (EP::EnumerativeProblem)(L::Loop)
    (S,P) = base_fibre(EP)
    loop = L.P
    for l in loop[1:end]
        (S,P) = (EP((S,P),l),l)
    end
    l = base_fibre(EP)[2]
    (S,P) = (EP((S,P),l),l)
    return(S)
end

"""
    permutation!(EP::EnumerativeProblem, ML::MonodromyLoop)

    Computes the permutation of the solutions of the enumerative problem `EP` induced by the monodromy loop `ML`.
"""
function permutation!(EP::EnumerativeProblem,ML::MonodromyLoop)
    if ML.sigma !== nothing
        return ML.sigma
    else
        sigma = permutation(EP, ML)
        ML.sigma = sigma
    end

end


function permutation(EP::EnumerativeProblem, ML::MonodromyLoop)
    g = numerical_function(solutions(base_fibre(ML)),EP(ML))
    if is_permutation(g,degree(EP)) == false
        @vprintln("The loop does not yield a valid permutation of the solutions of the enumerative problem.")
        @vprintln(" The function is: ",g)
        return nothing
    else
        return(perm(g))
    end 
end



#When sampling many monodromy elements, one can make use of the parallelization in HC.jl by choosing loops all of the form
#   a->b[i]->c->a where a is the base parameter, and c is another fixed parameter.
#   This requires 2*N+1 fibres to be tracked. 1 from a->c, and n_monodromy_loops from a->b[i] and c->b[i]
function large_monodromy_sample(F::System, bf::Fibre; n_monodromy_loops::Int = 50, monodromy_loop_scaling = 1.0, permutations_only = true)
    k = n_parameters(F)
    c = monodromy_loop_scaling.*randn(ComplexF64,k)
    b = [monodromy_loop_scaling.*randn(ComplexF64,k) for i in 1:n_monodromy_loops]
    a = parameters(bf)
    S1_bij = solve(F, bf[1]; start_parameters = bf[2], target_parameters=b)
    S2 = solve(F, bf[1]; start_parameters = bf[2], target_parameters = c)
    S2_bij = solve(F,S2; start_parameters = c, target_parameters = b)
    one_line_perms = [numerical_function(solutions(S2_bij[i][1]),solutions(S1_bij[i][1])) for i in 1:n_monodromy_loops]
    indices_of_valid_permutations = findall(x->is_permutation(x,length(bf[1])),one_line_perms)
    @vprintln("# Loops computed:            ",n_monodromy_loops)
    @vprintln("# Valid permutations:        ",length(indices_of_valid_permutations))
    @vprintln("# Unique valid permutations: ",length(unique(one_line_perms[indices_of_valid_permutations])))
    sampled_loops = [[a,bb,c,a] for bb in b]
    ML_bucket = Vector{MonodromyLoop}([])
    for i in indices_of_valid_permutations
        push!(ML_bucket,MonodromyLoop(bf,sampled_loops[i],perm(one_line_perms[i])))
    end
    return(ML_bucket)
end


const MONODROMY_SAMPLE = EnumerativeProperty{Vector{MonodromyLoop}}("monodromy_sample")

monodromy_sample(EP::EnumerativeProblem; kwargs...) = MONODROMY_SAMPLE(EP; kwargs...)


#TODO: When EnumerativeSolver is coded, that should be the input property
const large_monodromy_sample_datum = AlgorithmDatum(
    name = "large_monodromy_sample",
    description = "Sample of [n_monodromy_loops] random loops of scaling [monodromy_loop_scaling]",
    input_properties = [SYSTEM, BASE_FIBRE],
    default_kwargs = Dict{Symbol,Any}(:n_monodromy_loops=>50,:monodromy_loop_scaling=>1.0),
    output_property = MONODROMY_SAMPLE,
    reliability = :numerical_and_random
)

ALGORITHM_DATA[large_monodromy_sample]=large_monodromy_sample_datum



function subgroup(perms::Vector{PermGroupElem})
    H,_  = sub(symmetric_group(degree(perms[1])),perms)
    return(H)
end

"""
    group_generated_by_monodromy_loops(MLS::Vector{MonodromyLoop})

    Computes the subgroup of the symmetric group generated by the permutations induced by the monodromy loops in `MLS`.
"""
function group_generated_by_monodromy_loops(MLS::Vector{MonodromyLoop})
    mg = minimal_generating_set(subgroup(unique([ML.sigma for ML in MLS])))
    return(subgroup(mg))
end

const MONODROMY_GROUP = EnumerativeProperty{PermGroup}("monodromy group")


"""
    monodromy_group(EP::EnumerativeProblem; kwargs...)

Compute the monodromy group of the enumerative problem `EP`. Relevant keywords arguments, passed to `large_monodromy_sample`, are:
- `n_monodromy_loops`: The number of monodromy loops to sample (default: 50).
- `monodromy_loop_scaling`: The scaling factor for the monodromy loops (default: 1.0).
"""
monodromy_group(EP::EnumerativeProblem; kwargs...) = MONODROMY_GROUP(EP; kwargs...)
"""
    galois_group(EP::EnumerativeProblem; kwargs...)

Compute the Galois group of the enumerative problem `EP`. This is an alias for `monodromy_group`.
"""
galois_group(EP::EnumerativeProblem; kwargs...) = MONODROMY_GROUP(EP; kwargs...)




const group_generated_by_monodromy_loops_datum = AlgorithmDatum(
    name = "group_generated_by_monodromy_loops",
    description = "Computes the subgroup "*short(Dummit04)*" generated by monodromy permutations of the enumerative problem "*short(Jordan1870)*".",
    input_properties = [MONODROMY_SAMPLE],
    output_property = MONODROMY_GROUP,
    reliability = :certified,
    citations = [Dummit04, Jordan1870]
)

ALGORITHM_DATA[group_generated_by_monodromy_loops]=group_generated_by_monodromy_loops_datum


function is_decomposable(EP::EnumerativeProblem)
    G = galois_group(EP)
    if is_transitive(G)
        return !is_primitive(G)
    else
        return false
    end
end


function orbits(G::PermGroup)
    O = map(collect,oscar_orbits(G))
    return O
end
