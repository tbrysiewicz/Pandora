
function cr_num(z1,z2,z3,z4)
    (z1-z2)*(z3-z4)
 end

function cr_den(z1,z2,z3,z4)
    ((z1-z3)*(z2-z4))
 end


function distinct_entries(v)
	for i in 1:length(v)-1
		for j in i+1:length(v)
			if norm(v[i]-v[j])<0.000000001
				return(false)
			end
		end
	end
	return(true)
end


function CRsystem(B)
    @var Z[1:8]
    @var a[1:5]
    polys=[]
    for i in 1:5
        t = B[i]
        push!(polys,cr_num(Z[t]...)-cr_den(Z[t]...)*a[i])
    end
    push!(polys,Z[1]-10)
    push!(polys,Z[2]-1)
    push!(polys,Z[3]-5)
    start_pts = vcat([10.0,1.0,5.0],randn(ComplexF64,5))
    start_params = []
    for i in 1:5
        push!(start_params,cr_num(start_pts[B[i]]...)/cr_den(start_pts[B[i]]...))
    end
    
    F = System(polys,variables=Z,parameters=a)
    E = EnumerativeProblem(F)
    println("Bezout:",bezout(E))
    println("BKK:",bkk(E))
    P = randn(ComplexF64,5)
    S = HomotopyContinuation.solve(E.F,target_parameters=P,start_system = :total_degree, show_progress=false)
    S = Result(filter(x->x.return_code==:success,collect(S)))
    #Pandora.degree(E; Method = "Explicit")
    #S = base_solutions(E)
    fa=findall(x->distinct_entries(x),solutions(S))
    println("    only ",length(fa)," out of ",length(S)," are configurations of distinct points")
    if length(fa)>0
    E.BaseFibre = (Result(S[fa]),P)
    return(E)
    else
        "There are no solutions"
    return(nothing)
    end
    #MS = monodromy_solve(F,[start_pts],start_params)
    #return(MS)
end
