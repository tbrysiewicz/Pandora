export 
    MonodromyLoop,
    Path,
    monodromy,
    large_monodromy_sample,
    monodromy_group,
    galois_group,
    monodromy_sample

mutable struct MonodromyLoop
    F :: Fibre                      # Gives sols and ordering on them. (S,P)
    P :: Vector{Vector{ComplexF64}} # A list of parameters starting
    sigma :: Union{PermGroupElem,Nothing}
end


function Base.show(io::IO, ML::MonodromyLoop)
    print(io,"A monodromy ",length(ML.P),"-gon with permutation ",ML.sigma)
end


function MonodromyLoop(EP::EnumerativeProblem,P::Vector{Vector{ComplexF64}}) 
    F = base_fibre(EP)
    MonodromyLoop(F,P,nothing)
end


#p is interpretted as a function where p(i) = p[i]. But p(i) may be not injective or
#  surjective, in which case the following function should return false.
function is_valid_permutation(p::Union{Vector{Int},Vector{Union{Nothing,Int64}}},d::Int64) :: Bool
    u=unique(p) 
    return(!(length(u)<d || in(nothing,u))) #nothing indicates not surjective, length indicates not injective
end

function numerical_bijection(S1::Vector{T} where T, S2::Vector{T} where T; tol::Float64 = 1e-6) 
    one_line = [findfirst(x->isapprox(x,s),S1) for s in S2]
    return(one_line)
end

function check_and_fix_monodromy_input(fibre::Fibre, loop::Vector{Vector{T}})  where T
    if T != ComplexF64 #Make sure parameters are complex
        loop=Vector{Vector{ComplexF64}}(loop) #Parameter type: changing parameters in loop to complex floats
    end

    if isapprox(parameters(fibre),first(loop))==false #If input loop is doesn't start at fibre force it to
        pushfirst!(loop,parameters(fibre)) #Inconsistent loop: starting at given fibre parameter
    end

    if isapprox(first(loop),last(loop))==false #If input loop doesn't end at fibre, force it to
        push!(loop, first(loop)) #Incompelte loop: using start parameter as end parameter
    end
    return(loop)
end

function monodromy(EP::EnumerativeProblem, fibre::Fibre, loop::Vector{Vector{T}} where T)
    loop = check_and_fix_monodromy_input(fibre,loop)

    (S,P) = fibre
    for l in loop[2:end]
        (S,P) = (EP((S,P),l),l)
    end
    
    g = numerical_bijection(solutions(fibre),S)
    return is_valid_permutation(g,degree(EP)) ? perm(g) : g
end

monodromy(EP::EnumerativeProblem,loop::Vector{Vector{T}} where T) =monodromy(EP,base_fibre(EP),loop)
monodromy(EP::EnumerativeProblem,ML::MonodromyLoop) = monodromy(EP,ML.F,ML.P)
function monodromy!(EP::EnumerativeProblem,ML::MonodromyLoop)
    m = monodromy(EP,ML.F,ML.P)
    if typeof(m) <: Perm
        ML.sigma = m
    end
    return(m)
end

perm!(EP::EnumerativeProblem, ML::MonodromyLoop) = monodromy!(EP,ML)



#When sampling many monodromy elements, one can make use of the parallelization in HC.jl by choosing loops all of the form
#   a->b[i]->c->a where a is the base parameter, and c is another fixed parameter.
#   This requires 2*N+1 fibres to be tracked. 1 from a->c, and N from a->b[i] and c->b[i]
function large_monodromy_sample(S::System, F::Fibre; kwargs...)
    EP = EnumerativeProblem(S; populate = false)
    update_base_fibre!(EP,F)
    return(large_monodromy_sample(EP;kwargs...))
end

function large_monodromy_sample(EP::EnumerativeProblem; N::Int = 50, loop_scaling = 1.0, permutations_only = false)
    k = n_parameters(EP)
    c = loop_scaling.*randn(ComplexF64,k)
    b = [loop_scaling.*randn(ComplexF64,k) for i in 1:N]
    bf = base_fibre(EP)
    a = parameters(bf)
    S1_bij = solve(EP, base_fibre(EP), b)
    S2 = solve(EP,c)
    S2_fibre = (S2,c)
    S2_bij = solve(EP,S2_fibre,b)
    one_line_perms = [numerical_bijection(S2_bij[i],S1_bij[i]) for i in 1:N]
    indices_of_valid_permutations = findall(x->is_valid_permutation(x,degree(EP)),one_line_perms)
    println("# Loops computed:            ",N)
    println("# Valid permutations:        ",length(indices_of_valid_permutations))
    println("# Unique valid permutations: ",length(unique(one_line_perms[indices_of_valid_permutations])))
    if permutations_only
        return([perm(p) for p in one_line_perms[indices_of_valid_permutations]])
    else
        sampled_loops = [[a,bb,c,a] for bb in b]
        ML_bucket = Vector{MonodromyLoop}([])
        for i in indices_of_valid_permutations
            push!(ML_bucket,MonodromyLoop(bf,sampled_loops[i],perm(one_line_perms[i])))
        end
        return(ML_bucket)
    end
end


const MONODROMY_SAMPLE = EnumerativeProperty{Vector{MonodromyLoop}}("monodromy sample")

monodromy_sample(EP::EnumerativeProblem; kwargs...) = MONODROMY_SAMPLE(EP; kwargs...)


const large_monodromy_sample_datum = AlgorithmDatum(
    name = "Sample of [N] random monodromy loops of scaling [loop_scaling]",
    input_properties = [SYSTEM, BASE_FIBRE],
    default_kwargs = Dict{Symbol,Any}(:N=>50,:loop_scaling=>1.0),
    output_property = MONODROMY_SAMPLE,
    reliability = :numerical_and_random
)

ALGORITHM_DATA[large_monodromy_sample]=large_monodromy_sample_datum



function subgroup(perms::Vector{PermGroupElem})
    H,_  = sub(symmetric_group(degree(perms[1])),perms)
    return(H)
end

function group_generated_by_monodromy_loops(MLS::Vector{MonodromyLoop})
    return(subgroup([ML.sigma for ML in MLS]))
end

const MONODROMY_GROUP = EnumerativeProperty{PermGroup}("monodromy group")


"""
    monodromy_group(EP::EnumerativeProblem)

This function returns the `MONODROMY_GROUP` of an enumerative problem.

Note: Pandora.jl will automatically find an algorithm to compute the monodromy group of 
        EP. If you want to use a specific algorithm, call `algorithms_which_return(MONODROMY_GROUP)`
        If you have a preference of algorithm, run:
            `monodromy_group(EP::EnumerativeProblem; algorithm = )`

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
