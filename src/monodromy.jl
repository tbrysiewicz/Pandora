export 
    MonodromyLoop,
    Path,
    monodromy,
    large_monodromy_sample,
    monodromy_group,
    galois_group,
    monodromy_sample

mutable struct MonodromyLoop{T}
    F :: Fibre                      # Gives sols and ordering on them. (S,P)
    P :: Vector{Vector{T}}          # A list of parameters starting and ending at F[2]
    sigma :: Union{PermGroupElem,Nothing}
end

"""
    fibre(ML::MonodromyLoop)

    Returns the base fibre of the monodromy loop `ML`.
"""
function base_fibre(ML::MonodromyLoop)
    return ML.F
end

"""
    is_valid(ML::MonodromyLoop) 

    Checks if the monodromy loop `ML` is valid, meaning that it starts and ends at the base parameter of the fibre.
"""
function is_valid(ML::MonodromyLoop)
    base_param = fibre(ML)[2]   
    if isapprox(ML.P[1],base_param) && isapprox(ML.P[end],base_param) && base_param <: Vector{ComplexF64}
        return true
    else
        return false
    end
end

"""
    make_valid!(ML::MonodromyLoop) 

    Checks if the monodromy loop `ML` is valid, meaning that it starts and ends at the base parameter of the fibre.
"""
function make_valid!(ML::MonodromyLoop)
    if !is_valid(ML)
        base_param = fibre(ML)[2]
        if base_param <: Vector{Float64}
            ML.fibre = Fibre(ML.F.S, Vector{ComplexF64}(base_param)) # Convert to ComplexF64 
        end
        if !isapprox(ML.P[1],base_param)
            pushfirst!(ML.P, base_param) # Inconsistent loop: starting at given fibre parameter
        end
        if !isapprox(ML.P[end],base_param)
            push!(ML.P, base_param) # Incomplete loop: using start parameter as end parameter
        end
    end
    return(ML)
end

function Base.show(io::IO, ML::MonodromyLoop)
    if ML.sigma === nothing
        print(io,"A monodromy ",length(ML.P),"-gon. Use `monodromy!(EP,ML)` to compute the permutation.")
    else
        print(io,"A monodromy ",length(ML.P),"-gon with permutation ",ML.sigma)
    end
end


function MonodromyLoop(EP::EnumerativeProblem,P::Vector{Vector{ComplexF64}}) 
    F = base_fibre(EP)
    MonodromyLoop(F,P,nothing)
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


function monodromy!(EP::EnumerativeProblem,ML::MonodromyLoop)
    (S,P) = fibre(ML)
    for l in new_loop[2:end]
        (S,P) = (EP((S,P),l),l)
    end
    
    g = numerical_function(solutions(fibre),S)
    if is_permutation(g,degree(EP)) == false
        println("The loop does not yield a valid permutation of the solutions of the enumerative problem.")
        if typeof(loop) <: Vector{Vector{Float64}}
            println("  Are you sure you want to compute monodromy over the real numbers?")
        end
        return nothing
    else
        return(MonodromyLoop(fibre,new_loop,perm(g)))
    end 
end


function monodromy(EP::EnumerativeProblem, fibre::Fibre, loop::Vector{Vector{T}} where T)
    ML = MonodromyLoop(EP,fibre,loop)
    make_valid!(ML)
    return monodromy!(EP,ML)
end

monodromy(EP::EnumerativeProblem,loop::Vector{Vector{T}} where T) =monodromy(EP,base_fibre(EP),loop)

#=
function monodromy!(EP::EnumerativeProblem,ML::MonodromyLoop)
    m = monodromy(EP,ML.F,ML.P)
    if typeof(m) <: Perm
        ML.sigma = m
    end
    return(m)
end
=#

perm!(EP::EnumerativeProblem, ML::MonodromyLoop) = monodromy!(EP,ML)



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
    indices_of_valid_permutations = findall(x->is_permutation(x,length(bf)),one_line_perms)
    println("# Loops computed:            ",n_monodromy_loops)
    println("# Valid permutations:        ",length(indices_of_valid_permutations))
    println("# Unique valid permutations: ",length(unique(one_line_perms[indices_of_valid_permutations])))
    sampled_loops = [[a,bb,c,a] for bb in b]
    ML_bucket = Vector{MonodromyLoop}([])
    for i in indices_of_valid_permutations
        push!(ML_bucket,MonodromyLoop(bf,sampled_loops[i],perm(one_line_perms[i])))
    end
    return(ML_bucket)
end


const MONODROMY_SAMPLE = EnumerativeProperty{Vector{MonodromyLoop}}("monodromy sample")

monodromy_sample(EP::EnumerativeProblem; kwargs...) = MONODROMY_SAMPLE(EP; kwargs...)


#TODO: When EnumerativeSolver is coded, that should be the input property
const large_monodromy_sample_datum = AlgorithmDatum(
    name = "Sample of [n_monodromy_loops] random monodromy loops of scaling [monodromy_loop_scaling]",
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

function group_generated_by_monodromy_loops(MLS::Vector{MonodromyLoop})
    return(subgroup(unique([ML.sigma for ML in MLS])))
end

const MONODROMY_GROUP = EnumerativeProperty{PermGroup}("monodromy group")


"""
    monodromy_group(EP::EnumerativeProblem)

Compute the monodromy group of the enumerative problem `EP`.
"""
monodromy_group(EP::EnumerativeProblem; kwargs...) = MONODROMY_GROUP(EP; kwargs...)
galois_group(EP::EnumerativeProblem; kwargs...) = MONODROMY_GROUP(EP; kwargs...)




const group_generated_by_monodromy_loops_datum = AlgorithmDatum(
    name = "subgroup generated by monodromy sample",
    input_properties = [MONODROMY_SAMPLE],
    output_property = MONODROMY_GROUP,
    reliability = :heuristic
)

ALGORITHM_DATA[group_generated_by_monodromy_loops]=group_generated_by_monodromy_loops_datum


