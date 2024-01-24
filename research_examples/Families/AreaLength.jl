using Pandora
using HomotopyContinuation
using Oscar

function HeronFormula(a,b,c,A)
    (4*a*b-(a+b-c)^2-16*A)
end

function AreaLengthSystem(; positive_parameters=true)
    n = 5
    @var x[1:n,1:n]
    @var v[1:n,1:n,1:n]
    
    parameter_power=1
    if positive_parameters==true
        parameter_power=2
    end

    T = [[1,2,3],[1,2,4],[1,2,5],[1,3,4],
         [1,3,5],[1,4,5],[2,3,4],[2,3,5],
         [2,4,5],[3,4,5]]
    Eqs = [HeronFormula(
           x[i[1],i[2]],
           x[i[1],i[3]],
           x[i[2],i[3]],
           v[i[1],i[2],i[3]]^(parameter_power)) for i in T]
    F = System(Eqs, parameters = [v[i...] for i in T])
    E = EnumerativeProblem(F)
end


E = AreaLengthSystem()
d = degree(E)
G = galois_group(E)
B = minimal_block_reps(G)
order(G)
2^(32)*factorial(big(32))

function interval_sampler_sqrt(alpha)
    function sampler(n_samples,sample_length)
        #This uniformly samples from [alpha,1]
        P = [rand(sample_length).*(1-alpha).+alpha for p in 1:n_samples]
        #This takes square roots (inducing a non-uniform sample on [sqrt(alpha),1]) in preparation
        #  for squaring the parameters in the Area Length System
        P = [[sqrt(c) for c in p] for p in P]
    end
    return(sampler)
end

function interval_sampler(alpha)
    function sampler(n_samples,sample_length)
        #This uniformly samples from [alpha,1]
        P = [rand(sample_length).*(1-alpha).+alpha for p in 1:n_samples]
    end
    return(sampler)
end

function GramMatrix(x)
    [2*x[1] x[1]+x[2]-x[3] x[1]+x[4]-x[5] x[1]+x[7]-x[8]
    x[1]+x[2]-x[3] 2*x[2] x[2]+x[4]-x[6] x[2]+x[7]-x[9]
    x[1]+x[4]-x[5] x[2]+x[4]-x[6] 2*x[4] x[4]+x[7]-x[10]
    x[1]+x[7]-x[8] x[2]+x[7]-x[9] x[4]+x[7]-x[10] 2*x[7]]
end

function Euclidean(x)
    G = GramMatrix(x)
    println(G)
    e = LinearAlgebra.eigvals(G)
    println(e)
    for c in e
        if c<0
            println(c)
            return(false)
        end
    end
    return(true)
end


function Lorentzian(x)
    G = GramMatrix(x)
    println(G)
    e = LinearAlgebra.eigvals(G)
    pos = filter(x->x>0,e)
    if length(pos)==3
        return(true)
    end
    return(false)
end

function NEuclidean(E,R,P)
    length(filter(x->Euclidean(x),HomotopyContinuation.real_solutions(R)))
end
function NLorentzian(E,R,P)
    length(filter(x->Lorentzian(x),HomotopyContinuation.real_solutions(R)))
end

 v = [randn(Float64,4) for i in 1:5]
 EdgeOrder = [[1,2],[1,3],[2,3],[1,4],[2,4],[3,4],[1,5],[2,5],[3,5],[4,5]]
 X = [sum((v[e[1]]-v[e[2]]).^2) for e in EdgeOrder]
 G = GramMatrix(X)
 Euclidean(X)

Data = explore(E, [real_points_in_fibre, positive_points_in_fibre, NEuclidean, NLorentzian]; 
               sampler = interval_sampler(0.5), n_samples = 10000)

H = histogram([Data[2][1], Data[2][2]], labels=["Real" "Positive"], bins=0:1:64)
H = histogram([filter(x->x!=0,Data[2][4]),filter(x->x!=0,Data[2][3])], labels=["Real" "Positive"],bins=1:1:32)

