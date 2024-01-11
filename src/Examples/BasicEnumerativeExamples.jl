export
    univariate_family, #lower case if input is taken #camel case otherwise
    AreaLengthSystem





@doc raw"""
    univariate_family(d::Int64)

 Returns the enumerative problem of solving a generic univariate polynomial of degree $d$.
 # Examples
 ```jldoctest
 julia> E = univariate_family(5);

 julia> degree(E)
 Populating a base fibre of the enumerative problem
 5
 ```
 """
function univariate_family(d::Int64)
    @var x
    @var a[1:d+1]
    F = System([a[d+1]+sum([a[i]*x^i for i in 1:d])],variables=[x],parameters=vec(a))
    E = EnumerativeProblem(F)
    return(E)
end




function AreaLengthSystem()
    function HeronFormula(a,b,c,A)
        (4*a*b-(a+b-c)^2)-16*A
    end
    @var x[1:5,1:5]
    @var v[1:5,1:5,1:5]

    T = combinations(1:5,3)
    Eqs =[HeronFormula(x[i[1],i[2]],x[i[1],i[3]],x[i[2],i[3]],v[i[1],i[2],i[3]]) for i in T]
    F = System(Eqs,parameters = [v[i[1],i[2],i[3]] for i in T])
    E = EnumerativeProblem(F)
end




function PosSampler(n,k;low=0,high=1)
    c = high-low
    [[(rand(Float64)*c+low)^2 for j in 1:k] for i in 1:n]
end

function ALsampler1(n)
    PosSampler(n,10;low=0,high=1)
end
function ALsampler2(n)
    PosSampler(n,10;low=0.3,high=1)
end

function ALsampler3(n)
    PosSampler(n,10;low=0.5,high=1)
end

function ALsampler4(n) #randomly pick points to make a 4 simplex, and then pull the
                        #areas via Heron
    P = [randn(Float64,2) for i in 1:5]
    T = collect(combinations(P,3))
    E = [[norm(t[1]-t[2]),norm(t[1]-t[3]),norm(t[2]-t[3])] for t in T]
    println(E)
    function AreaFromEdges(a,b,c)
        (1/4)*sqrt((4*a^2*b^2-(a^2+b^2-c^2)^2))
    end
    A = [AreaFromEdges(e...) for e in E]
end
