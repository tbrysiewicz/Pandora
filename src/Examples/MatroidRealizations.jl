export 
	matrix_to_nonbases,
	matroid_realization_space,
	best_realizable_matroid




###################################################################################################
###########################CONSTRUCTING MATROID REALIZATION SPACES ################################
###################################################################################################

#Function: matroid_space_eqns
#Input: n = number of elements in a matroid, nonbases = nonbases of a rank 3 matroid
#Output: A tuple containing the equations of describing a matroid space and a general matrix from 
#        which the equations are derived from. 
function matroid_space_eqns(n::Int64, nonbases::Vector{Vector{Int64}}; rank = 3, reduced = false) 
	if reduced == true
		M = matroid_from_nonbases(nonbases,n)
		matroid_space_eqns(M; reduced = true)
	end
	@var x[1:rank-1,1:n]
	matrix = vcat([1 for i in 1:n]', x) 
	eqns = [det(matrix[1:(rank),c]) for c in nonbases] #Computing equations 
	return (System(eqns, variables = vcat(x...)), matrix)	
end
function matroid_space_eqns(M::Matroid; reduced = false)
	if reduced == true
		println("NOT FINISHED")
		RS = realization_space(M)
		HC_vars =  [Variable(x) for x in RS.ambient_ring.S]
		for v in RS.ambient_ring.S
			@var v
			push!(HC_vars,v)
		end
		return(HC_vars)
	end
	NB = Vector{Vector{Int}}(nonbases(M))
	n = length(M.groundset)
	return(matroid_space_eqns(n, NB; rank=rank(M)))
end 



@doc raw"""
    matroid_realization_space(n::Int,nonbases::Vector{Vector{Int}}) 
    matroid_realization_space(Matr::Matroid)
    
Returns the variety of the realization space of a rank three matroid over the complex numbers. If the matroid is not realizable over the complex numbers, nothing will be returned. 
# Examples
 ```jldoctest
 julia> matroid_realization_space(7, [[1,2,3],[1,4,5],[1,6,7]])
✓ Decomposing 1 witness sets    Time: 0:00:01
  Current status:                    
  Degrees of components of dim. 11:  8
There are 1 many component(s) of the realization space.
They are all 11-dimensional          1
Assigning the corresponding witness (superset) to the witness slot of the variety
A variety in C^14 of degree 8 and dimension 11.

julia> matroid_realization_space(fano_matroid())
✓ Decomposing 2 witness sets    Time: 0:00:01
  Current status:                   
  Degrees of components of dim. 9:  6
  Degrees of components of dim. 7:  2,2,2,2,2,2,2,2,2,2,2,2,2,2
Matroid is not realizable

```
"""
function matroid_realization_space(n::Int,nonbases::Vector{Vector{Int}}; reduced = false)  
	(eqns,matrix_format) = matroid_space_eqns(n, nonbases; reduced = reduced)
	if length(eqns) == 0
		println("ERROR: decide how to deal with unirational realization spaces")
	end
	matroid_variety = Variety(eqns)
	NID = nid(matroid_variety)
	sorted_nonbases = sort([sort(m) for m in nonbases])
	WS = witness_sets(NID)

	CorrectComponents = []

	for i in keys(WS) #Scroll through each dimension

		number_components_dim_i = length(WS[i]) #number of components with dimension i 
		
		for j in 1:number_components_dim_i #Scroll through each of the components of that dimension

			solution = sample(WS[i][j]) #sample a point on that component
			matrix_solution = HomotopyContinuation.evaluate(matrix_format, variables(matrix_format)=> solution) 
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
		return(matroid_variety,matrix_format)
	else
		println("ERROR: This code needs to be written - when correct components have different dimensions")
	end
end

function matroid_realization_space(Matr::Matroid; reduced = false)
	NB = Vector{Vector{Int}}(nonbases(Matr))
	n = length(Matr)
	return(matroid_realization_space(n, NB; reduced = reduced))
end 

###################################################################################################
###########################          HELPER FUNCTIONS              ################################
###################################################################################################

 #TODO: (for matrix_to_nonbases) Eventually switch to a reliable function from /AlgebraicMatroids
 @doc raw"""
 	matrix_to_nonbases(matrix::Union{Matrix{Int64},Matrix{Float64},Matrix{ComplexF64}}; tol = 10000) 
 
Returns the corresponding nonbases of $matrix$. 
# Examples
```jldoctest
julia> M = [1 0 1 4; 0 1 1 5; 0 0 0 6]
3×4 Matrix{Int64}:
1  0  1  4
0  1  1  5
0  0  0  6

julia> matrix_to_nonbases(M)
1-element Vector{Any}:
[1, 2, 3]

```
"""
function matrix_to_nonbases(matrix::Union{Matrix{Int64},Matrix{Float64},Matrix{ComplexF64}}; tol = 10000) 
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
 

###################################################################################################
###########################    VISUAL OPTIMIZATION FUNCTIONS       ################################
###################################################################################################

#this uses the matrix_format of a matroid realization space to create a function which scores a
#  matroid fibre's visual appeal by judging each matroids appeal in the fibre, and taking min
function create_matroid_fibre_visual_appeal(matrix_format)
	function mfva(F ::Tuple{Result,Vector{Float64}})
		real_sols= HomotopyContinuation.real_solutions(F[1]) #real solutions in fibre
       
		if real_sols != []
			N = length(real_sols) #N of real solutions 
			record_max = maximum([matroid_visual_appeal(real_sols[i],matrix_format) for i in 1:N]) #Take maximum score
		else 
				   println("there are no real solutions") 
				return(-Inf)
		end
		return(record_max)
	end
	return(mfva)
end

#matroid_visual_appeal: scores a matroid based on its visual appeal   
#  WARNING: this assumes the first row is all ones (see second line)
# Reshapes the solution S into a matrix M
function matroid_visual_appeal(S,matrix_format)
	M_projective = HomotopyContinuation.evaluate(matrix_format, variables(matrix_format)=> S) 		
	M = M_projective[2:3,:]
	line_intersection_points = line_intersection_points_of_matroid(M)
	new_M = hcat(M, line_intersection_points...) #Matrix with matroid points and line intersection points
	n = size(new_M,2) #the number of elements in our matroid
	point_pairs = collect(combinations(1:n,2)) #all possible combinations of pairs of points
	record_min = Inf 
	for j in 1:(length(point_pairs)) #scrolls through all possible point pairs 	
		p = new_M[:,(point_pairs[j][1])] #first point in the pair
		q = new_M[:,(point_pairs[j][2])] #second point in the pair
		x = norm(p-q) #norm between the 2 points in pair j 
		s =  -(1/1.851)*(x-0.1)*(x-2*sqrt(2)) #scoring that norm
		
		if s < record_min
			record_min = s #recording the minimum norm 
		end
	end   
	return(record_min)    
end

#This function is left in for me to compare how old score function without intersection points compares to the new one. 
function old_matroid_visual_appeal(S)		
	M = reshape(S, 2, convert(Int, length(S)/2)) #reshapes solution vector s into a matrix
	n = size(M,2) #the number of elements in our matroid
	point_pairs = collect(combinations(1:n,2)) #all possible combinations of pairs of points
	record_min = Inf 
	for j in 1:(length(point_pairs)) #scrolls through all possible point pairs 	
		p = M[:,(point_pairs[j][1])] #first point in the pair
		q = M[:,(point_pairs[j][2])] #second point in the pair
		x = norm(p-q) #norm between the 2 points in pair j 
		s =  -(1/1.851)*(x-0.1)*(x-2*sqrt(2)) #scoring that norm
		
		if s < record_min
			record_min = s #recording the minimum norm 
		end
	end   
	return(record_min)    
end

################################################################################################
############################### VISUAL APPEAL HELPER FUNCTIONS #################################
################################################################################################



#TODO: allow for other hyperplanes at infinity
#TODO: allow for points at infinity too!
function projective_to_affine(M::Matrix; tol = 0.00001) 
	(m,n) = size(M)
	for i in 1:n
		if abs(M[1,i])<tol
			println("This matrix represents projective configuration with points at infinity")
			return(nothing)
		end
		for j in 1:m
			M[j,i]=M[j,i]/M[1,i]
		end
	end #After this step, the matrix has first row equal to the all ones vector
	return(M[2:3,:])
end


#Finds the intersection point of 2 lines in the plane, where line 1 is given by points p1 and p2, and line 2 is given by points p3 and p4. 
function intersection_point_of_lines(p1, p2, p3, p4)
	a = (p2[2] - p1[2]) / (p2[1] - p1[1]) #slope of the first line given by p1 and p2
	b = (p4[2] - p3[2]) / (p4[1] - p3[1]) #slope of the second line given by p3 and p4
	x = (1 / (a - b)) * (p3[2] - p1[2] + (a * p1[1]) - (b * p3[1])) #x coordinate of the intersection point
	y =  (a*(x-p1[1])) + p1[2] #y coordinate of the intersection point 
	intersection_point = vcat(x,y)
	return(intersection_point)
end 

#Outputs a vector of vectors where each vector contains all points that lie on a single line in the matroid 
function matroid_lines_nonbases(nonbases) 
	L = [] #collection of lines describing the matroid  
	
	for i in 1:length(nonbases) 
		nb = nonbases[i]
		l = nb #vector that has points that lie on a line with the points in nb (initializing with nb)
		nonbases_without_nb = filter(x -> x != nb, nonbases)
	
		for j in 1:length(nonbases_without_nb)
			common_elements = intersect(nb, nonbases_without_nb[j])
			
			if length(common_elements) == 2 
				l = sort(unique(vcat(l, nonbases_without_nb[j]))) #adding any elements not already in l but in nonbases_without_nb and sorting it 
			end 
		end
		
		if !(l in L) && !(isempty(l)) #adding any unique lines to the collection of lines L
			push!(L,l)
		end 
	end 
	return(L)	
end 

#Finds all intersection points of the lines of a matroid drawing in the plane, where intersection points that already corresond to points of the matroid are not included in output. 
function line_intersection_points_of_matroid(M)
	n = size(M,2) #the number of elements in our matroid
	projective_M = vcat([1 for i in 1:n]', M) 
	nb = matrix_to_nonbases(projective_M) #nonbases of the matrix 
	lines_as_nonbases = matroid_lines_nonbases(nb) #lines of matroid described by nonbases 
	k = length(lines_as_nonbases)
	line_pairs = collect(combinations(1:k,2)) #all possible line pairs 
	intersection_points = []

	for i in 1:length(line_pairs)
		first_nb_in_line_pair = lines_as_nonbases[line_pairs[i][1]]
		second_nb_in_line_pair = lines_as_nonbases[line_pairs[i][2]]
		common_elements = intersect(first_nb_in_line_pair, second_nb_in_line_pair) #Checking if a point in matroid is an intersection of line pairs 

		if isempty(common_elements) 
			p1 = M[:,first_nb_in_line_pair[1]] #column vector of M corresponding to the first entry in the first nonbases in the line pair 
			p2 = M[:,first_nb_in_line_pair[2]]
			p3 = M[:,second_nb_in_line_pair[1]]
			p4 = M[:,second_nb_in_line_pair[2]]
			intersection_point = intersection_point_of_lines(p1,p2,p3,p4)
			push!(intersection_points, intersection_point)
			
		end 
	end 
	return(intersection_points)
end 





################################################################################################
############################### WRAPPER FOR OPTIMIZATION PROBLEM################################
################################################################################################


@doc raw"""
    best_realizable_matroid(n:: Int64, nonbases:: Vector{Vector{Int64}}; n_trials::Int64 = 50, n_samples::Int64 = 10)
    best_realizable_matroid(M::Matroid; n_trials = 50, n_samples=10)
    
Returns a real matrix repersentative of a rank three matroid. This matrix repersentative has the best visual appeal when drawing the matroid in the plane. If the matroid is not realizable or if there are no real matrix repersentatives, then nothing is returned.  

Matrices are sampled from the realization space of the matroid, and then optimized based on a visual appeal score. n_trials corresponds to the number of trials the optmizer will take and n_samples corresponds to the number of samples the optmizer will take on each trial. 
# Examples
 ```jldoctest
julia> best_realizable_matroid(non_fano_matroid(), n_trials = 2)
✓ Decomposing 2 witness sets    Time: 0:00:00
  Current status:                   
  Degrees of components of dim. 9:  6
  Degrees of components of dim. 8:  3,3,22,3,1
There are 1 many component(s) of the realization space.
They are all 8-dimensional
Assigning the corresponding witness (superset) to the witness slot of the variety
-----------------------------------------------Trial: 1
Number of samples: 10
Number of successful tracks: 10
there are no real solutions
there are no real solutions
there are no real solutions
there are no real solutions
there are no real solutions
there are no real solutions
there are no real solutions
there are no real solutions
Current Score: -0.6008297846505287
     Status: Improved Record Score
     Improvement: 15.952363197997954
Number of fibres computed:10
Number of non-taboo fibres:10
-----------------------------------------------Trial: 2
Number of samples: 10
Number of successful tracks: 10
there are no real solutions
there are no real solutions
there are no real solutions
there are no real solutions
there are no real solutions
there are no real solutions
there are no real solutions
Current Score: 0.08141457545739453
     Status: Improved Record Score
     Improvement: 0.6822443601079232
Number of fibres computed:10
Number of non-taboo fibres:10
2×7 Matrix{Float64}:
 -0.71116   0.162516  0.420361  -0.691239  -0.700836  0.439747  0.937398
  0.669725  1.11295   1.24375    0.23219    0.44298   1.39895   1.71426

```
"""
function best_realizable_matroid(n:: Int64, nonbases:: Vector{Vector{Int64}}; n_trials::Int64 = 50, n_samples::Int64 = 10, shotgun_hill_climb::Int64 = 3, verbose = true)

	(V,matrix_format) = matroid_realization_space(n, nonbases) #variety associated to matroid space equations
	matroid_fibre_visual_appeal = create_matroid_fibre_visual_appeal(matrix_format)
	if !isnothing(V) #if the matroid has a realization space
		E = EnumerativeProblem(V) #Turning V into an ennumerative problem 
		path_solutions = HomotopyContinuation.solutions(E.BaseFibre[1]) #all solutions from each pathresult of the base fibre
		best_path_solution = path_solutions[argmax([matroid_visual_appeal(s,matrix_format) for s in path_solutions])] #finding the solution (matroid) in path_solutions that has the highest score
		P = base_parameters(E)
		NewR = solve(system(E), best_path_solution; start_parameters = P, target_parameters =P) #converting best_path_solution to a result 
		E.BaseFibre = (NewR, P) #Changing base fibre of E so that it only tracks the path associated to best_path_solution

		best_score = -inf #Intializng best score for shot-gun hill climb
		S =  [] #initalizing S 

		for i in 1:shotgun_hill_climb
			OE = initialize_optimizer(E, matroid_fibre_visual_appeal)
			OE = optimize!(OE;n_trials=n_trials, n_samples=n_samples, verbose=verbose) #finding the set of solutions with the best configuration of points 
			current_score = OE.record_score 
			i = i +1 
			if current_score > best_score 
				best_score = current_score #updating the best score to be the current score 
				S = HomotopyContinuation.real_solutions(OE.record_fibre[1]) #taking the real solutions from the results 
			end 
		end 

		if !is_empty(S) #if the optimizer returns a real solution
			best_matroid_as_matrix =HomotopyContinuation.evaluate(matrix_format, variables(matrix_format)=> S[1])
			return(projective_to_affine(best_matroid_as_matrix))
		else 
			println("During optimization, we found no real matrix repersentatives.")
			return(nothing)
		end 

	else 
		println("The realization space is empty")
		return(nothing)
	end 
end 

function best_realizable_matroid(M::Matroid; n_trials = 50, n_samples=10)
	NB = Vector{Vector{Int}}(nonbases(M))
	n = length(M.groundset)
	best_realizable_matroid(n,NB;n_trials=n_trials,n_samples=n_samples)
end
