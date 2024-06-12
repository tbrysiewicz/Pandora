export 
	matroid_space_eqns,
	points_to_matrix,
	matrix_to_nonbases,
	find_real_points,
	random_real_linear_space, 
	matroid_realization_space


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
function matroid_space_eqns(n::Int, nonbases::Vector{Vector{Int}}) 
	@var x[1:2,1:n] #variables are (x_{1,i},x_{2,i}) for each point
	matrix = vcat([1 for i in 1:n]', x) 
	eqns = [det(matrix[1:3,c]) for c in nonbases] #Computing equations #TCB: TODO: try a non-determinantal set of equations at some point. 
	return(System(eqns,variables=vcat(x...)))	
end 


function matroid_space_eqns(Matr::Matroid)
	NB = Vector{Vector{Int}}(nonbases(Matr))
	n = length(Matr)
	return(matroid_space_eqns(n, NB))
end 


#Function : solution_to_matrix 
#Input: points = vector of point entries  
#		d = dimension of space points lie in
#Output: m = matrix where columns are the points 
#TODO: generalize this function. What if the solution corresponds to a differently shaped matrix?
function points_to_matrix(points::Union{Vector{Float64},Vector{ComplexF64}}, d::Int) 
	n = div(length(points),d) #Number of columns 
 	M = reshape(points, d, n)  #Matrix of points
	return(M)
 end 

#Function : matrix_to_nonbases 
#Input: matrix = matirx   
#Output: nonbases = Non-Bases of matrix
#TODO: Eventually switch to a reliable function from /AlgebraicMatroids
function matrix_to_nonbases(matrix::Union{Matrix{Float64},Matrix{ComplexF64}}; tol = 10000) 
	nonbases= []
	r = rank(matrix,0.0001)
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
	matroid_variety = Variety(matroid_space_eqns(n, nonbases))
	NID = nid(matroid_variety)
	sorted_nonbases = sort([sort(m) for m in nonbases])
	WS = witness_sets(NID)

	CorrectComponents = []

	for i in keys(WS) #Scroll through each dimension

		number_components_dim_i = length(WS[i]) #number of components with dimension i 
		
		for j in 1:number_components_dim_i #Scroll through each of the components of that dimension

			solution = sample(WS[i][j]) #sample a point on that component
			matrix_solution = vcat([1 for i in 1:n]', points_to_matrix(solution, 2)) #turn that point into a matrix (projective)
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
function find_real_points(V::Variety,N::Int; limit = 100)
	D = ambient_dimension(V) 
	d = D -Pandora.dim(V) #dimension of linear subspace is codim of variety
	W = witness_set(V) 
	realpts = [] #initalizing list of real points 
	i = 0 #initalizing counter

	while i < limit 

		L = random_real_linear_space(d, D) #random linear subspace
		newW = witness_set(W,L) #witness set of W using a random linear subspace
		realpts = vcat(realpts, HomotopyContinuation.real_solutions(results(newW))) #list of real points off of NewW
		length(realpts)>= N && break #break if we have found at least N real points 
		i = i+1
	end 
	Nrealpts = vcat(realpts[1:N]) #collecting N real points off of witness set
	return(Nrealpts)
end 

#Function : random_real_linear_space
#Input: d = dimension of linear space, D = dimension of ambient space 
#Output: L = random affine linear space of dimension d
function random_real_linear_space(d::Int,D::Int)
	L = rand_subspace(D, dim=d, real = true, affine = true)
	return(L)
end
     

#Function : matrix_repersentative 
#Input: n = number of points, d = dimension of space, nobases = nonbases of matroid  
#Output: M = matrix repersentative of matroid 
function matrix_repersentative(n::Int, d::Int, nonbases::Vector{Vector{Int}})
	V = matroid_realization_space(n,nonbases)
	if V == nothing
		println("There is no matrix repersentative for this matroid.") 
		return(nothing)
	
	else
		pt= find_real_points(V,1)
		M = points_to_matrix(pt[1],d)
		return(M)
	end 
end



#These are some score functions I wrote, I'm still working on assigning a specific score to the matroid, instead of just computing data about it 
#(e.g. slopes or distances). I will put it all under one master score function later. 
function distance_score(M::Matrix{Float64})
	n = size(M,2) #number of points (columns of matrix)
	point_pairs = collect(combinations(1:n, 2)) #all point pair combinations
	d = []

	for j in 1:length(point_pairs) 	
		u = M[:,(point_pairs[j][1])]
		v = M[:,(point_pairs[j][2])]
		push!(d, distance(u, v, EuclideanNorm())) #putting the distance between all pairs of points into a vector
	end 						
	return(maximum(abs,d), minimum(abs,d)) #returning the largest and smallest distance between a pair of points 
end 

function slope_score(M::Matrix{Float64}, nonbases::Vector{Vector{Int}})
	x = M[1,:] 
	y = M[2,:]
	for i in 1:size(M,2) 
		point = i 
		NB_with_i = filter(x -> point in x, nonbases) #list of all nonbases with point i in it 
		slopes = [] #initializing vector of slopes of lines all going through point i 

		for j in 1:length(NB_with_i)
			slopes = push!(slopes, (y[NB_with_i[j][2]] - y[NB_with_i[j][1]]) / (x[NB_with_i[j][2]] - x[NB_with_i[j][1]])) #add slopes of each line to list
			line_pairs = collect(combinations(1:length(NB_with_i))) 

			for k in length(line_pairs)
				u = slopes[line_pairs[k][1]]
				v = slopes[line_pairs[k][2]]
				d = u-v 
			end 

		end 

	end

end 