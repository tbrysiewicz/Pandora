export 
	matrix_to_nonbases,
	matroid_realization_space,
	best_realizable_matroid


#TODO List (of functions)
#matroid_space_eqns
#	- Should I generalize this so it works for matroids of different rank? 
#
#points_to_matrix 
#
#matroid_realization_space
#	- Should I generalize this so it works for matroids of different rank?
#	- How to deal with the matroid given by nonbases = []. As the code errors out since the NID cannot be computed. 
#	- Discuss the issue of a trace test failure when the NID is computed with Taylor. (Motivated by the matroid that Taylor emailed to me)
#	- Discuss what Taylor is referring to in the last block of code that has println("ERROR: This code needs to be written")
#	
#find_N_real_points 
#	- Maybe delete?? No longer used in my code.
#
#matroid_fibre_visual_appeal
#	- Now that I have changed best_realizable_matroid so that the base fibre of E so only tracks the path associated to best_path_solution, so the fibre only has one matroid solution, should I delete this?? I would have to change the code of best_realizable_matroid to use matroid_visual_appeal as a score. 
#	
#matroid_visual_appeal
#	- Add type for the input 
#	- If I delete matroid_fibre_visual_appeal, I would need to change this up, specifcally add a message if the solution isn't real. 
#	- This needs to be edited for including the scoring of intersection points 
#
#best_realizable_matroid 
#
#RECENTLY ADDED FUNCTIONS:
#intersection_point_of_lines
#
#matroid_lines_nonbases
#
#line_intersection_points_of_matroid
#	
#GENERAL COMMENTS 
#- For the matroid Taylor emailed me, it was sent as M = [[3,8,9],[4,7,8],[1,2,6,8],[8,11,12],[5,6,10],[9,10,11],[2,4,5,9],[1,9,12],[1,3,5,7],[3,10,12]]. There was a nonbases of 4 elements, should I make a function that breaks this up, and then call it in my other functions?? As currently, my functions won't take this input as the nonbases are not the same length.
#


#Function: matroid_space_eqns
#Input: n = number of elements in a matroid, nonbases = nonbases of a rank 3 matroid
#Output: A tuple containing the equations of describing a matroid space and a general matrix from which the equations are derived from. 
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
#Output: m = matrix where columns are the points. If projective = true, a row of ones is added to the top of the matrix. 
function points_to_matrix(points::Union{Vector{Float64},Vector{ComplexF64}}, d::Int; projective=true) 
	n = div(length(points),d) #Number of columns 
 	M = reshape(points, d, n)  #Matrix of points
	if projective == true
		M = vcat([1 for i in 1:n]', M)
	end
	return(M)
 end 


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
function matroid_fibre_visual_appeal(F :: Tuple{Result,Vector{Float64}}) #scoring the variety 
	real_sols= HomotopyContinuation.real_solutions(F[1]) #real solutions in fibre
       
	if real_sols != []
		N = length(real_sols) #N of real solutions 
		record_max = maximum([matroid_visual_appeal(real_sols[i]) for i in 1:N]) #Take maximum score
	else 
       		println("there are no real solutions") 
			return(-Inf)
	end
	return(record_max)
end
       

#sc: scores a matroid solution -- what if they choose a different initialized basis than [0,1], [-1,0], [0,1]  
# Reshapes the solution S into a matrix M. Then scrolls through all possible point pairs and calculated the
#norm between each pair of points and scores that norm on the function s. The minimum score across all point pairs is the
#score of our matrix
function matroid_visual_appeal(S)		
	M = reshape(S, 2, convert(Int, length(S)/2)) #reshapes solution vector s into a matrix
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
function best_realizable_matroid(n:: Int64, nonbases:: Vector{Vector{Int64}}; n_trials::Int64 = 50, n_samples::Int64 = 10, shot_gun_hill_climb::Int64 = 3)

	V = matroid_realization_space(n, nonbases) #variety associated to matroid space equations
	if !isnothing(V) #if the matroid has a realization space
		E = EnumerativeProblem(V) #Turning V into an ennumerative problem 
		path_solutions = HomotopyContinuation.solutions(E.BaseFibre[1]) #all solutions from each pathresult of the base fibre
		best_path_solution = path_solutions[argmax([matroid_visual_appeal(s) for s in path_solutions])] #finding the solution (matroid) in path_solutions that has the highest score
		P = base_parameters(E)
		NewR = solve(system(E), best_path_solution; start_parameters = P, target_parameters =P) #converting best_path_solution to a result 
		E.BaseFibre = (NewR, P) #Changing base fibre of E so that it only tracks the path associated to best_path_solution

		best_score = -inf #Intializng best score for shot-gun hill climb
		S =  [] #initalizing S 

		for i in 1:shot_gun_hill_climb
			OE = initialize_optimizer(E, Pandora.matroid_fibre_visual_appeal)
			OE = optimize!(OE;n_trials=n_trials, n_samples=n_samples) #finding the set of solutions with the best configuration of points 
			current_score = OE.record_score 
			i = i +1 
			if current_score > best_score 
				best_score = current_score #updating the best score to be the current score 
				S = HomotopyContinuation.real_solutions(OE.record_fibre[1]) #taking the real solutions from the results 
			end 
		end 

		if !is_empty(S) #if the optimizer returns a real solution
			best_matroid_as_matrix = points_to_matrix(S[1],2;projective=false)
			return(best_matroid_as_matrix)
		else 
			println("There are no real matrix repersentatives.")
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




