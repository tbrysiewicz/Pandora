

function secant_type(S)
    (t1,t2,s1,s2) = S
    t = [t1,t2]
    s = [s1,s2]
    #first check if t[1] and t[2] are approximately real 
    conj_pairs = 0
    real_pairs = 0 
    if isapprox(t[1], conj(t[2]), atol=1e-6)
        conj_pairs += 1
    #else if t[1] and t[2] are approximately real
    elseif norm(imag(t[1])) < 1e-6 && norm(imag(t[2])) < 1e-6
        real_pairs += 1
    end
    #now check if s[1] and s[2] are approximately real
    if isapprox(s[1], conj(s[2]), atol=1e-6)
        conj_pairs += 1
    elseif norm(imag(s[1])) < 1e-6 && norm(imag(s[2])) < 1e-6
        real_pairs += 1
    end
    if conj_pairs == 2
        return [0,0,1]
    elseif real_pairs == 2
        return [1,0,0]
    elseif conj_pairs == 1 && real_pairs == 1
        return [0,1,0]
    else
        return [0,0,0]
    end
end

function total_secant_type(S)
    if length(S)==0
        return [-1,-1,-1]
    end
    sum(map(secant_type,S)).//4
end

function confirm_type(EP, P)
    testp = vcat(eachcol(reshape(P, 4, 4)')...)
    S = EP(testp)
    return Pandora.total_secant_type(S)
end

function load_types()
    found_types = []
    associated_fibres = []
    #parse the odd lines of the file "DoubleSecants.txt" to get the types and the even lines to get the associated fibres
    open(pwd()*"/OutputFiles/DoubleSecants.txt", "r") do io
        lines = filter(!isempty, strip.(readlines(io)))
        if isodd(length(lines))
            throw(ArgumentError("Malformed DoubleSecants.txt: expected pairs of lines (type, fibre)."))
        end
        for i in 1:2:length(lines)
            type_line = lines[i]
            if occursin("[", type_line)
                type = Core.eval(@__MODULE__, Meta.parse(type_line))
            else
                type = parse.(Int, split(type_line, ","))
            end
            push!(found_types, type)

            fibre_line = lines[i+1]
            if occursin("(", fibre_line)
                fibre = Core.eval(@__MODULE__, Meta.parse(fibre_line))
            else
                fibre = parse.(Float64, split(fibre_line, ","))
            end
            push!(associated_fibres, fibre)
        end
    end
    return found_types, associated_fibres
end




#=

using LinearAlgebra


function secant_type(S)
    (t1,t2,s1,s2) = S
    t = [t1,t2]
    s = [s1,s2]
    #first check if t[1] and t[2] are approximately real 
    conj_pairs = 0
    real_pairs = 0 
    if isapprox(t[1], conj(t[2]), atol=1e-6)
        conj_pairs += 1
    #else if t[1] and t[2] are approximately real
    elseif norm(imag(t[1])) < 1e-6 && norm(imag(t[2])) < 1e-6
        real_pairs += 1
    end
    #now check if s[1] and s[2] are approximately real
    if isapprox(s[1], conj(s[2]), atol=1e-6)
        conj_pairs += 1
    elseif norm(imag(s[1])) < 1e-6 && norm(imag(s[2])) < 1e-6
        real_pairs += 1
    end
    if conj_pairs == 2
        return [0,0,1]
    elseif real_pairs == 2
        return [1,0,0]
    elseif conj_pairs == 1 && real_pairs == 1
        return [0,1,0]
    else
        return [0,0,0]
    end
end

function total_secant_type(S)
    if length(S)==0
        return [-1,-1,-1]
    end
    sum(map(secant_type,S))
end

function confirm_type(EP, P)
    testp = vcat(eachcol(reshape(P, 4, 4)')...)
    S = EP(testp)
    return Pandora.total_secant_type(S)
end

function load_types()
    found_types = []
    associated_fibres = []
    #parse the odd lines of the file "DoubleSecants.txt" to get the types and the even lines to get the associated fibres
    open(pwd()*"/OutputFiles/DoubleSecants.txt", "r") do io
        lines = filter(!isempty, strip.(readlines(io)))
        if isodd(length(lines))
            throw(ArgumentError("Malformed DoubleSecants.txt: expected pairs of lines (type, fibre)."))
        end
        for i in 1:2:length(lines)
            type_line = lines[i]
            if occursin("[", type_line)
                type = Core.eval(@__MODULE__, Meta.parse(type_line))
            else
                type = parse.(Int, split(type_line, ","))
            end
            push!(found_types, type)

            fibre_line = lines[i+1]
            if occursin("(", fibre_line)
                fibre = Core.eval(@__MODULE__, Meta.parse(fibre_line))
            else
                fibre = parse.(Float64, split(fibre_line, ","))
            end
            push!(associated_fibres, fibre)
        end
    end
    return found_types, associated_fibres
end
(found_types, associated_fibres) = load_types()

EP = twisted_cubic_bisecants() #Below, we use the 10-to-1 map instead of the 40-to-1 map. 
S = base_fibre(EP)[1]
P = base_fibre(EP)[2]
S = sort(S, lt=(x,y)->real(sum(x))<real(sum(y)))
S = S[[1,5,9,13,17,21,25,29,33,37]]
EP.knowledge=EP.knowledge[[1,2]]
Pandora.know!(EP,Pandora.BASE_FIBRE,(S,P))



for i in 1:10000
	P1 = [randn(Float64,16) for i in 1:1000]
    P2 = [[rand() < 0.33 ? Float64(1.0) : (rand() < 0.20/0.67 ? Float64(0.0) : randn()*rand(1:4)) for _ in 1:16] for _ in 1:1000]
    P3 = [[rand() < 0.33 ? Float64(1.0) : (rand() < 0.20/0.67 ? Float64(0.0) : randn()*rand(-4:4)) for _ in 1:16] for _ in 1:1000]
    P4 = [[rand() < 0.33 ? Float64(-1.0) : (rand() < 0.20/0.67 ? Float64(1.0+rand()*0.001) : Float64(0)) for _ in 1:16] for _ in 1:1000]
	P5 = [[randn(Float64) for j in 1:16] for i in 1:1000]
    P6 = [vcat([1.0,0.0,0.0,0.0],randn(Float64,12).*rand(1:6)) for i in 1:10000]
    P7 = [vcat(eachcol(reshape(p,4,4)')...) for p in P6] 
    P8 = [vcat([1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],randn(Float64,8).*rand(1:6)) for i in 1:10000]
    P9 = [vcat(eachcol(reshape(p,4,4)')...) for p in P8] 
    #scale the 1st, 6th, 11th and 16th entries of P7 and P9 by 10^(rand(-3:3)) to get more variety in the types we find
    P10 = [vcat(eachcol(reshape(p,4,4)'*reshape(p,4,4))...) for p in P7]
    for i in eachindex(P7); P7[i] = P7[i] .* [2.0^rand(-4:4) for _ in 1:16]; P7[1][1] = 1.0; end
    for i in eachindex(P7); P9[i] = map(abs,P7[i]); P9[1][1] = 1.0; end
    #Make sure P7 and P9 are lists of 16-dimensional vectors and not Vector{Base.Generator...}
    #Now set P to be the union of P1 to P4
    P = vcat(P1, P2, P3, P4, P5,P7,P9,P10)
    SP = EP(P)
	FIB = [Fibre(solutions(sp[1]),sp[2]) for sp in SP]
	TST = [total_secant_type(f[1]) for f in FIB]
	U = unique(TST)
	for u in U
		if (u ∉ found_types)
			uindex = findfirst(tst -> tst==u,TST)
			push!(associated_fibres,FIB[uindex])
			push!(found_types,u)
			open(pwd()*"/OutputFiles/DoubleSecants.txt", "a") do io
    				println(io, string(u))
    				println(io, string(FIB[uindex]))
			end
			println("Found new type!!!                                     ",u)
		end
	end
	
	 for t in sort(filter(x->sum(x)==10,found_types))
                   println(t)
               end
	println("Total number of types found:",length(found_types))
	println("Total number of types in last layer:",length(filter(x->sum(x)==10,found_types)))
end
=#

function Segre_CC(n)

    C = Pandora.cube(n,0,1)
    L=Vector{Vector{Int}}(Pandora.lattice_points(C))
    
    @var x[1:n]
    
    function mon(v,e)
        prod([v[i]^e[i] for i in 1:length(v)])
    end
    
    @var c[1:n+1, 1:2^n]
    
    @var lambda
    
    Eqs = [sum([c[1,j]*mon(x,L[j]) for j in 1:2^n])-lambda]
    println(Eqs[1])
    for i in 1:n
        f = sum([c[i+1,j]*mon(x,L[j]) for j in 1:2^n])-lambda*x[i]
        println(f)
        push!(Eqs,f)
    end
        
    F = System(Eqs,variables=vcat(vec(x),lambda),parameters=vec(c))
    E = EnumerativeProblem(F)

end



function Segre_CC_Decomposed(n)

    C = Pandora.cube(n,0,1)
    L=Vector{Vector{Int}}(Pandora.lattice_points(C))
    
    @var x[1:n]
    
    @var a[1:n,1:n]
    
    @var lambda
    
    Eqs = [prod([x[i]-2^i for i in 1:n])-lambda]
    println(Eqs[1])
    for i in 1:n
        f = prod([x[j]-a[i,j] for j in 1:n])-lambda*x[i]
        println(f)
        push!(Eqs,f)
    end
        
    F = System(Eqs,variables=vcat(vec(x),lambda),parameters=vec(a))
    E = EnumerativeProblem(F)

end

#-335.3553125015627, -29.546992252608817, 424.716804286756, 395.7250779531599, 163.8221461732173, -376.367090236099, 114.51013200021491, 123.71361996223806, 121.43004013473022, 249.2370481503624, -664.2962691218697, 91.47394525909013, 125.26568380283216, 178.43609867104203, -47.02164234479464, 326.0222889909182, -47.25811693197321, 178.4248006018958, -11.838772090897306, 15.04646200317575, 137.32222476842432, -390.7907535502272, 138.84250423494643, -0.6262558419921923, 177.72879145684857, 19.201420675934283, 1.7984122353115104, -457.09178750046, -273.7468376099906, -147.17373390612602]


O = initialize_real_optimizer(E)
O.sampler.n_samples=1000
improve!(O)
improve!(O)
improve!(O)
O.sampler.n_samples=60
for i in 1:10
    O.optimizer_data.step=0
    O.optimizer_data.steps_no_progress=0
    O.optimizer_data.steps_no_major_progress=0
    optimize!(O)
end