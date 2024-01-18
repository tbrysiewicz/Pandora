
using LinearAlgebra
using HomotopyContinuation

function GenericEDMatrix(n)
    @var x[1:n+1,1:n+1]
    M = [x[min(i,j),max(i,j)]*(i!=j) for i in 1:n+1, j in 1:n+1]
end

###########################################################Basic Constructors
#Takes as input a Euclidean Distance Matrix
# D = (d_ij^2) where d_ij is the distance between
# point p_i and p_j; equivalently, the squared
# edge lengths of any n-simplex with vertices (p_i)
function CayleyMengerDeterminant(D;flag="nothing")
    n=size(D)[1]
    newRow = [1 for i in 1:n]
    newCol = vcat([1 for i in 1:n],0)
    M = hcat(vcat(D,newRow'),newCol)
    C = (-1)^(n)//((factorial(n-1))^2*2^(n-1))
    if flag=="NoConstant"
    	return(det(M))
    else
    	return(C*det(M))
    end
end


function CMDeterminants(D,C; flag = "nothing")
    n = size(D)[1]
    #C = collect(filter(x->length(x)>1,collect(combinations(1:n))))
    DetBucket=[]
    for c in C
        push!(DetBucket,CayleyMengerDeterminant(D[c,c];flag=flag))
    end
    return(DetBucket)
end


function ConstructHeronSystemGeneral(B,n)
      X = GenericEDMatrix(n)

      Vals = CMDeterminants(X,B)
      Vars = Vector{Variable}(unique(vcat(X...))[2:end])

      @var y[1:length(B)]
      F = System([Vals[i]-y[i] for i in 1:length(B)],variables = Vars,parameters = y)

      M = GenericEDMatrix(n)
      SVars = Vector{Variable}(unique(vcat(M...))[2:end])
      Msubbed = subs(M,SVars=>randn(ComplexF64,length(SVars)))
      StartSol = Vector{ComplexF64}(unique(vcat(Msubbed...))[2:end])
      StartParams = CMDeterminants(Msubbed,B)
      return(F,StartSol,StartParams)
end

CandidateBases = [[[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 3], [1, 4], [2, 3], [2, 4]],
[[1, 2], [1, 2, 3], [1, 3], [1, 4], [2, 4], [3, 4]],
[[1, 2], [1, 2, 3, 4], [1, 3], [1, 4], [2, 3], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [1, 4], [2, 3]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [2, 3], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [1, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 3], [1, 4], [2, 3, 4], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [2, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 3], [2, 3, 4], [2, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [1, 3, 4], [2, 3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3], [1, 4], [2, 3]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3], [1, 4], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 4], [2, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3], [2, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3, 4], [2, 3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [1, 3, 4], [1, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [1, 4], [2, 3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [1, 3, 4], [2, 3]],
[[1, 2], [1, 2, 3], [1, 3, 4], [1, 4], [2, 3, 4], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [1, 3, 4], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [2, 3, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3], [1, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3], [1, 4], [2, 3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3], [2, 3]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3, 4], [1, 4], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3], [2, 3, 4], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3, 4], [2, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3], [1, 3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3], [2, 3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3, 4], [1, 4], [2, 3, 4]]]

BI = [1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1]

HeronProblems = []
for i in 1:35
    if BI[i]==1
        E = EnumerativeProblem(ConstructHeronSystemGeneral(CandidateBases[i],3)[1])
        degree(E)
        push!(HeronProblems,E)
    else
        push!(HeronProblems,nothing)
    end
end

GaloisGroups = [] 
for i in 1:35
    if BI[i]==1
        E = HeronProblems[i]
        push!(GaloisGroups,galois_group(E;nloops = 20))
    else
        push!(GaloisGroups,nothing)
    end
end

function EDfrom3Sol(s)
    M = [0 s[1] s[2] s[3]; s[1] 0 s[4] s[5]; s[2] s[4] 0 s[6]; s[3] s[5] s[6] 0]
end

allvols = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4],[1,2,3],[1,2,4],[1,3,4],[2,3,4],[1,2,3,4]]

function VolVectorFromSols(s)
    D = EDfrom3Sol(s)
    CM = CMDeterminants(D,allvols; flag = "nothing")
    return(CM)
end

CSG = []
for i in 1:35
    if BI[i]==1
        E = HeronProblems[i]
        push!(CSG,coordinate_symmetry_group(E;F = VolVectorFromSols))
    else
        push!(CSG,nothing)
    end
end

for i in 1:35
    if BI[i]==1
        if GaloisGroups[i]!=CSG[i]
            println(i,"    ",describe(GaloisGroups[i]),"     ",describe(CSG[i]))
        end
    end
end




























BASES = CandidateBases[findall(x->x==1,BI)]
HERONPROBLEMS = [EnumerativeProblem(ConstructHeronSystemGeneral(B,3)[1]) for B in BASES]
[Degree(E;Retry=true,Method="Explicit") for E in HERONPROBLEMS]

E = HERONPROBLEMS[7]


function SolToCayleyMengerMatrix(s)
    M = [0 1 1 1 1; 1 0 s[1] s[2] s[3]; 1 s[1] 0 s[4] s[5]; 1 s[2] s[4] 0 s[6]; 1 s[3] s[5] s[6] 0]
end

S = (solve_over_params(E,[randn(Float64,6) for i in 1:100000]))

function triangleInequality(T)
    for a in T
        if a<0
            return(false)
        end
    end
    t = [sqrt(T[i]) for i in 1:3]
    if t[1]+t[2]>t[3] && t[1]+t[3]>t[2] && t[2]+t[3]>t[1]
        return(true)
    else
        return(false)
    end
end
function SimplexSolutions(s)
    RS = real_solutions(s)
    RealizeSimplices = []
    for p in RS
        Triangles = [[p[1],p[2],p[4]],[p[1],p[3],p[5]],[p[2],p[3],p[6]],[p[4],p[5],p[6]]]
        SatisfiesTriangleInequalities=true
        for T in Triangles
            if triangleInequality(T)==false
                SatisfiesTriangleInequalities=false
            end
        end
#        println("Triangle Inequalities: ",SatisfiesTriangleInequalities)
        M = SolToCayleyMengerMatrix(p)
#        println(M)
#        println(det(M))
        if SatisfiesTriangleInequalities==true && det(M)>0
            push!(RealizeSimplices,p)
            println("Realizes a simplex!")
        else
#            println(".")
        end

    end
    return(RealizeSimplices)
end

function CountSimplexSolutions(E,S,P)
    return(CountSimplexSolutions(S))
end

function CountSimplexSolutions(s)
    SS = SimplexSolutions(s)
    return(length(SS))
end


AllData = []
for i in 1:length(BASES)
    E = HERONPROBLEMS[i]
    Data = Explore(E,[RealPointsInFibre,CountSimplexSolutions,PositivePointsInFibre];n_samples=10000)
    push!(AllData,Data)
end

function DataHistogram(i; Condition=false)
    B = BASES[i]
    DATA = AllData[i]
    E = HERONPROBLEMS[i]
    NewData = DATA
    if Condition==true
        #Indicator = findall(x->x!=0,DATA[2][1])
        #NewData = [DATA[1],[d[Indicator] for d in DATA[2]]]
        NewData = [DATA[1], [filter(x->x!=0,d) for d in DATA[2]]]
    end
    DATA=NewData
    H1 = plot(labels=DATA[1][1])
    H2 = plot(labels=DATA[1][3])
    H3 = plot(labels=DATA[1][2])
    if length(DATA[2][1])>0
        H1=histogram(DATA[2][1],labels=DATA[1][1],color=:blue)
    end
    if length(DATA[2][2])>0
        H3=histogram(DATA[2][2],labels=DATA[1][2],color=:red)
    end
    if length(DATA[2][3])>0
        H2=histogram(DATA[2][3],labels=DATA[1][3],color=:green)
    end
    H = [H1,H2,H3]
    P=plot(title=string(B), grid = false, showaxis = false, bottom_margin = -5Plots.px)
    Q = plot(P,H..., layout = @layout([A{0.01h}; [B C D]]),size=(1000,300))
    return(Q)
end
