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
    Eqs =[HeronFormula(x[i[1],i[2]],x[i[1],i[3]],x[i[2],i[3]],v[i[1],i[2],i[3]]^2) for i in T]
    push!(Eqs,v[1,2,3]-1)
    push!(Eqs,v[1,2,4]-1)
    F = System(Eqs,parameters = [v[i[1],i[2],i[3]] for i in filter(x->x!=[1,2,3] && x!=[1,2,4]  ,collect(T))])
    E = EnumerativeProblem(F)
end


function TangentCirclesToCubics()
    @var s, t, r # Circle centered at (s, t) with radius r
    @var u[1:3], v[1:3] # Three points of tangency (u[1], v[1]), ..., (u[3], v[3])
    @var a[1:10], b[1:10], c[1:10] # Three fixed cubics with coefficients defined by a[1], ..., a[10]; ... ; c[1], ..., c[10]
    
    #=
    Set up the polynomial system defining tritangent circles
    =#
    
    f_1 = (u[1] - s)^2 + (v[1] - t)^2 - r # the point (u[1], v[1]) lies on the circle
    
    f_2 = (u[2] - s)^2 + (v[2] - t)^2 - r # the point (u[2], v[2]) lies on the circle
    
    f_3 = (u[3] - s)^2 + (v[3] - t)^2 - r # the point (u[3], v[3]) lies on the circle
    
    f_4 = a[1]*u[1]^3  + a[2]*u[1]^2*v[1] + a[3]*u[1]*v[1]^2 + a[4]*v[1]^3 + a[5]*u[1]^2 + a[6]*u[1]*v[1] + a[7]*v[1]^2 + a[8]*u[1] + a[9]*v[1] + a[10] 
    
    f_5 = b[1]*u[2]^3  + b[2]*u[2]^2*v[2] + b[3]*u[2]*v[2]^2 + b[4]*v[2]^3 + b[5]*u[2]^2 + b[6]*u[2]*v[2] + b[7]*v[2]^2 + b[8]*u[2] + b[9]*v[2] + b[10] 
    
    f_6 = c[1]*u[3]^3  + c[2]*u[3]^2*v[3] + c[3]*u[3]*v[3]^2 + c[4]*v[3]^3 + c[5]*u[3]^2 + c[6]*u[3]*v[3] + c[7]*v[3]^2 + c[8]*u[3] + c[9]*v[3] + c[10] 
    
    f_7 = det([differentiate(f_1, [u[1], v[1]]) differentiate(f_4, [u[1], v[1]])]) # Circle and conic are tangent at (u[1], v[1])
    
    f_8 = det([differentiate(f_2, [u[2], v[2]]) differentiate(f_5, [u[2], v[2]])]) # Circle and conic are tangent at (u[2], v[2])
    
    f_9 = det([differentiate(f_3, [u[3], v[3]]) differentiate(f_6, [u[3], v[3]])]) # Circle and conic are tangent at (u[3], v[3])
    
    # Polynomial system in variables  [u[1], v[1], u[2], v[2], u[3], v[3], s, t, r] and parameters [a[1], ..., a[6], b[1], ..., b[6], c[1], ..., c[6]]
    
    F = System([f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9], variables = [u[1], v[1], u[2], v[2], u[3], v[3], s, t, r], 
    parameters = [a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10],  
                  b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10],
                  c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10]]);
    
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
