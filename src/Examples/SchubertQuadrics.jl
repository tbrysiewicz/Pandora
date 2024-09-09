export
    SchubertQuadrics
function SchubertQuadrics(α,β,γ)
    Schuberts_Triangle = [[1,3,9,17,21,21,17,9,3,1],
                           [2,6,18,34,42,24,18,6,2],
                            [4,12,36,68,68,36,12,4],
                              [8,24,72,104,72,24,8],
                              [16,48,112,112,48,16],
                                  [32,80,128,80,32],
                                    [56,104,104,56],
                                        [80,104,80],
                                            [92,92],
                                               [92]]
    if α+β+γ != 9
        return("α+β+γ must equal 9")
    end
    enumerative_count=Schuberts_Triangle[β+1][α+1]
    @var p[1:4,1:α] l[1:6,1:β] h[1:4,1:γ]
    @var D x[1:4,1:4]
    #We represent quadrics by 4x4 symmetric matrices
    X = Matrix{Expression}(Symmetric(x))
    #Below we construct the second and third exterior power of X
    cols₂ =[[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
    Λ₂X = [det(X[I,J]) for I in cols₂, J in cols₂]

    cols₃ = [[2,3,4],[1,3,4],[1,2,4],[1,2,3]]
    S=[[1,0,0,0] [0,-1,0,0] [0,0,1,0] [0,0,0,-1]]
    Λ₃X= S*[det(X[I,J]) for I in cols₃, J in cols₃]*S
    
    Point_Conditions=[p'*X*p for p in eachcol(p)]
    Line_Conditions=[l'*Λ₂X*l for l in eachcol(l)]
    Plane_Conditions=[h'*Λ₃X*h for h in eachcol(h)]
    
    Affine_Chart=sum(randn(Float64,10).*unique(vec(X)))-1
    Det_Value=det(X)-D
    
    params=vcat(vec(p),vec(l),vec(h))
    Equations=System(vcat(Point_Conditions,
                          Line_Conditions,
                          Plane_Conditions,
                          Affine_Chart,
                          Det_Value), 
                          parameters=params)
    E=EnumerativeProblem(Equations)
    return(E)
end