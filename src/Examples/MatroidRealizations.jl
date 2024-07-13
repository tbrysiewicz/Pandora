export 
	matroid_space_eqns,
	points_to_matrix,
	matrix_to_nonbases,
	find_real_points,
	matroid_realization_space,
	Matroid_score,
	best_realizable_matroid


#TODO List
#-Take an oscar matroid as input 		
#		Alex: this is done
#-Put types in the inputs 
#		Alex: this is done 
#-Run tests (make another file with queries whose answers you know)
#		Alex: I would like to run a few more tests to really make sure it's solid 
#-Make sure the names of functions are exactly what they do 								
#            (e.g. solution to matrix is a general function name which does something
#                  very specific: it makes a 2xn matrix.)
#-Can you upgrade this to handle matroids of points in higher-dimensional spaces (other than the plane?)
#		Alex: I renamed solution_to_matrix and I did update this function, so it takes an 	
#			  input d::the dimension of space the points are in. I'm thinking I should change matroid_space_eqns and matroid_realization_space
#			  to also take points that are not just in the plane. 
#
#-Write a function which is called													
#   function find_real_points(V::Variety,N::Int; limit = 100)
#  which will compute (at most) "limit" - many witness sets over real linear spaces
#  and collect the real points along the way until N real points have been found. 
#-You will need to figure out a way to make a "random" linear space. Make that into a function
#  too: function random_linear_space(d::Int,D::Int) # make random linear space of dimension d in an ambient space of dimension D
#  use the HomotopyContinuation.jl type "LinearSpace"  				
#	Alex: This function is written. HomotopyContinuation had a bulit in function to produce a random linear space but I still wrote the random
#		  linearspace function using it. I'm not sure if I should just delete this since it only uses HomotopyContinuation's built in function?
#
# -Continue writing the function which scores the 'visual appeal' of a 2xn matrix (along w/ dependencies of interest) 
#          whose cols correspond to points
#	Alex: I'm working away at this. 
#
#- Write a function which draw the point and lines
# function draw_matroid_representative(M::Matrix{Float64},nonbases=Vector{Vector{Int}})
#  use scatter and use plot. Put it in the visualization folder in a file called
#   MatroidVisualization.jl
#	Alex: this is done
#
# -After writting this, you can write something like
# function draw_matroid_representative(M::Matrix{Float64},Matr::Matroid)
#    NB = nonbases(Matr)
#    return(draw_matroid_representative(M,NB))
# end
#	Alex: this is done
#
#
# -Together, by the end of the week it would  be great to have a function which samples 1000
#   real configurations representing a matroid, and picks the one that "looks" the best and
#   plots it. 
# -It would be nice if all of the relevant functions which one could reasonably expect
#  to call on an object of type "Oscar.Matroid" can be called in that way. 
# -If you really want to push hard this week, take a look at Optimization/CoreObjects at the "Score"
#   constructions there. 



#Function: matroid_space_eqns  
#Input: n = number of points,
#       NonBases = points we want collinear/ nonbases of matroid  
#Output: Eqns = determinant equations defining a matroid 
#		 m = matrix repersenting the realization space, with 3 columns forming a basis set as triangle (1,1,0),(1,-1,0) and (1,0,1) 

#function matroid_space_eqns(n::Int, nonbases::Vector{Vector{Int}}) 
	#bases = filter(x -> !(x in nonbases), collect(combinations(1:n,3))) #finding the bases
	#k = filter(x ->!(x in bases[1]), 1:n) #finding all points not in the first basis
	#j = 1 #initalizing counter
	#m = zeros(3,n) #initalizing matrix of all zeros 
	#X = [] #initalizing list of variables 
		
	#for i in 1:n 
		#if i in k 
			#@var x[1:2,i]
			#X = push!(X,x)
			#m = hcat(m[:, 1:(i-1)], vcat(1,x), m[:, (i+1):end])
		#end 
			
		#if i in bases[1]
			#M = [1 1 1;-1 1 0; 0 0 1]
			#m = hcat(m[:, 1:(i-1)], M[:,j], m[:, (i+1):end])
			#j = j+1 
		#end 
	#end
	
	#eqns = [det(m[1:3,c]) for c in nonbases] #Computing equations #TCB: TODO: try a non-determinantal set of equations at some point. 
	#return(System(eqns, variables=vcat(X...)), m)
#end 

function matroid_space_eqns(n::Int64, nonbases::Vector{Vector{Int64}}) 
	@var x[1:2,1:n]
	matrix = vcat([1 for i in 1:n]', x) 
	eqns = [det(matrix[1:3,c]) for c in nonbases] #Computing equations 
	return (System(eqns, variables = vcat(x...)), matrix)	
end
function matroid_space_eqns(M::Matroid)
	NB = Vector{Vector{Int}}(nonbases(M))
	n = length(M.groundset)
	return(matroid_space_eqns(n, NB))
end 

#Function : solution_to_matrix 
#Input: points = vector of point entries  
#		d = dimension of space points lie in
#Output: m = matrix where columns are the points 
#TODO: generalize this function. What if the solution corresponds to a differently shaped matrix?
function points_to_matrix(points::Union{Vector{Float64},Vector{ComplexF64}}, d::Int; projective=true) 
	n = div(length(points),d) #Number of columns 
 	M = reshape(points, d, n)  #Matrix of points
	if projective == true
		M = vcat([1 for i in 1:n]', M)
	end
	return(M)
 end 

#Function : matrix_to_nonbases 
#Input: matrix = matirx   
#Output: nonbases = Non-Bases of matrix
#TODO: Eventually switch to a reliable function from /AlgebraicMatroids
function matrix_to_nonbases(matrix::Union{Matrix{Float64},Matrix{ComplexF64}}; tol = 10000) 
	nonbases= []
	r = rank(matrix,1.0/tol)
	n = size(matrix, 2) #number of columns
	C = collect(combinations(1:n, r)) #collection of all possible bases 
	
	for i in 1:length(C)
		c=C[i]
		submatrix = matrix[:, c] #All possible submatrices
		submatrix_conditioning = cond(submatrix) 
	
		if abs(submatrix_conditioning) > 10000
			push!(nonbases, c) 
		end	
	end
	return nonbases	
end 
 

#Function : matroid_realization_space
#Input: n= number of points, nonbases = points we want collinear/ nonbases of a matrix  
#Output: matroidwitnessresults = Witness set of a matroid
function matroid_realization_space(n::Int,nonbases::Vector{Vector{Int}})  
	eqns = matroid_space_eqns(n, nonbases)
	matroid_variety = Variety(eqns[1])
	NID = nid(matroid_variety)
	sorted_nonbases = sort([sort(m) for m in nonbases])
	WS = witness_sets(NID)

	CorrectComponents = []

	for i in keys(WS) #Scroll through each dimension

		number_components_dim_i = length(WS[i]) #number of components with dimension i 
		
		for j in 1:number_components_dim_i #Scroll through each of the components of that dimension

			solution = sample(WS[i][j]) #sample a point on that component
			matrix_solution = HomotopyContinuation.evaluate(eqns[2], variables(eqns[2])=> solution) 
			nonbases_component = matrix_to_nonbases(matrix_solution) #compute the corresponding matroid
			sorted_nonbases_component= sort([sort(m) for m in nonbases_component])
                 
			if sorted_nonbases_component == sorted_nonbases #Adding witness set component of dimension i to vector if non-bases are equal to input 
				push!(CorrectComponents,[i,j])
			end 	
		end 
		
	end 

	if length(CorrectComponents)==0
		println("Matroid is not realizable")
		return(nothing)
	end

	if length(unique([x[1] for x in CorrectComponents])) == 1 #If the correct components all belong to the same dimension
		println("There are ",length(CorrectComponents)," many component(s) of the realization space.")
		println("They are all ", CorrectComponents[1][1],"-dimensional")
		println("Assigning the corresponding witness (superset) to the witness slot of the variety")
		assign_component!(matroid_variety,CorrectComponents[1][1],[cc[2] for cc in CorrectComponents])
		return(matroid_variety)
	else
		println("ERROR: This code needs to be written")
	end
end


function matroid_realization_space(Matr::Matroid)
	NB = Vector{Vector{Int}}(nonbases(Matr))
	n = length(Matr)
	return(matroid_realization_space(n, NB))
end 

#Function : find_real_points
#Input: V = variety, N = number of real points to find   
#Output: Nrealpts = list containing N real points from V
function find_N_real_points(V::Variety,N::Int; limit = 100)
	D = ambient_dimension(V) 
	d = D -Pandora.dim(V) #dimension of linear subspace is codim of variety
	W = witness_set(V) 
	realpts = [] #initalizing list of real points 
	i = 0 #initalizing counter

	while i < limit 

		L = rand_subspace(D, dim=d, real = true, affine = true)
		newW = witness_set(W,L) #witness set of W using a random linear subspace
		realpts = vcat(realpts, HomotopyContinuation.real_solutions(results(newW))) #list of real points off of NewW
		length(realpts)>= N && break #break if we have found at least N real points 
		i = i+1
	end 
	Nrealpts = vcat(realpts[1:N]) #collecting N real points off of witness set
	return(Nrealpts)
end 


     
#matroid_fibre_visual_appeal : scores a set of matroid solutions (scores a fiber) 
#input: F = ()
function matroid_fibre_visual_appeal(F :: Tuple{Result,Vector{Float64}}) #scoring the variety 
	real_sols= HomotopyContinuation.real_solutions(F[1]) #real solutions in fibre
       
	if real_sols != []
		N = length(real_sols) #N of real solutions 
		record_max = maximum([matroid_visual_appeal(real_sols[i]) for i in 1:N]) #Take maximum score
	else 
       		println("there are no real solutions") #will this happen? And what should I put?
			return(-Inf)
	end
	return(record_max)
end
       

#sc: scores a matroid solution -- what if they choose a different initialized basis than [0,1], [-1,0], [0,1]  
# Reshapes the solution S into a matrix M. Then scrolls through all possible point pairs and calculated the
#norm between each pair of points and scores that norm on the function s. The minimum score across all point pairs is the
#score of our matrix
function matroid_visual_appeal(S)		
	m = reshape(S, 2, convert(Int, length(S)/2)) #reshapes solution vector s into a matrix
	M = hcat(m, [0,1], [-1,0], [0,-1]) #adds set bases vectors to m 
	n = size(M,2) #the number of elements in our matroid
	point_pairs = collect(combinations(1:n,2)) #all possible combinations of pairs of points
	record_min = Inf 
	for j in 1:(length(point_pairs)) #scrolls through all possible point pairs 	
		p = M[:,(point_pairs[j][1])] #first point in the pair
		q = M[:,(point_pairs[j][2])] #second point in the pair
		x = norm(p-q) #norm between the 2 points in pair j 
		s =  -(1/1.851)*(x-0.1)*(x-2*sqrt(2)) #scoring that distance (I have no idea how big a marker is)
		
		if s < record_min
			record_min = s #recording the minimum distance 
		end
	end   
	return(record_min)    
end


#Function : best_realizable_matroid
#Input: n = number of points, nonbases = nonbases of matroid, N = number of steps optimization function will take   
#Output: A planar graph of the matroid with the best point configuration for visual appeal 

function best_realizable_matroid(n:: Int64, nonbases:: Vector{Vector{Int64}}; n_trials::Int64 = 50, n_samples::Int64 = 10)

	V = matroid_realization_space(n, nonbases) #variety associated to matroid space equations
	if !isnothing(V) #if the matroid has a realization space
		E = EnumerativeProblem(V) #Turning V into an ennumerative problem 
		OE = initialize_optimizer(E,matroid_fibre_visual_appeal)
		OE = optimize!(OE;n_trials=n_trials, n_samples=n_samples) #finding the set of solutions with the best configuration of points 
		S = HomotopyContinuation.real_solutions(OE.record_fibre[1]) #taking the real solutions from the results 
		if !is_empty(S) #if the optimizer returns a real solution
			best_matroid_in_fibre = S[argmax([matroid_visual_appeal(s) for s in S])] #vector of  scores of all solutions in S
			best_matroid_as_matrix = points_to_matrix(best_matroid_in_fibre,2;projective=false)
			P = draw_matroid_representative(best_matroid_as_matrix, nonbases)
			return(P)
		else 
			println("There are no real solutions")
			return(nothing)
		end 

	else 
		return(nothing)
	end 
end 

function best_realizable_matroid(M::Matroid; n_trials = 50, n_samples=10)
	NB = Vector{Vector{Int}}(nonbases(M))
	n = length(M.groundset)
	best_realizable_matroid(n,NB;n_trials=n_trials,n_samples=n_samples)
end



