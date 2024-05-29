export 
	matroid_space_eqns,
	matroid_realization_space

	

#TODO List
#-Take an oscar matroid as input
#-Put types in the inputs
#-Run tests (make another file with queries whose answers you know)
#-Make sure the names of functions are exactly what they do 
#            (e.g. solution to matrix is a general function name which does something
#                  very specific: it makes a 2xn matrix.)
#-Can you upgrade this to handle matroids of points in higher-dimensional spaces (other than the plane?)
#-Write a function which is called
#   function find_real_points(V::Variety,N::Int; limit = 100)
#  which will compute (at most) "limit" - many witness sets over real linear spaces
#  and collect the real points along the way until N real points have been found. 
# You will need to figure out a way to make a "random" linear space. Make that into a function
#  too: function random_linear_space(d::Int,D::Int) # make random linear space of dimension d in an ambient space of dimension D
#  use the HomotopyContinuation.jl type "LinearSpace" 
# Continue writing the function which scores the 'visual appeal' of a 2xn matrix (along w/ dependencies of interest) 
#          whose cols correspond to points
# Write a function which draw the point and lines
# function draw_matroid_representative(M::Matrix{Float64},nonbases=Vector{Vector{Int}})
#  use scatter and use plot. Put it in the visualization folder in a file called
#   MatroidVisualization.jl
# After writting this, you can write something like
# function draw_matroid_representative(M::Matrix{Float64},Matr::Matroid)
#    NB = nonbases(Matr)
#    return(draw_matroid_representative(M,NB))
# end
# Together, by the end of the week it would  be great to have a function which samples 1000
#   real configurations representing a matroid, and picks the one that "looks" the best and
#   plots it. 
# It would be nice if all of the relevant functions which one could reasonably expect
#  to call on an object of type "Oscar.Matroid" can be called in that way. 
# If you really want to push hard this week, take a look at Optimization/CoreObjects at the "Score"
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


#Function : solution_to_matrix 
#Input: soln = solution vector  
#Output: m = matrix of a row of ones and solution
#TODO: generalize this function. What if the solution corresponds to a differently shaped matrix?
function solution_to_matrix(soln::Union{Vector{Float64},Vector{ComplexF64}}) 
	n = div(length(soln),2) #Number of columns 
 	matrix_without_ones = reshape(soln, 2, n)  #Matrix of points without the row of ones 
 	M = vcat([1 for i in 1:n]' , matrix_without_ones) #Adding the row of 1s to make M
	return(M)
 end 



#Function : matrix_to_nonbases 
#Input: matrix = matirx   
#Output: nonbases = Non-Bases of matrix
#TODO: Eventually switch to a reliable function from /AlgebraicMatroids
function matrix_to_nonbases(matrix; tol = 10000) 
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
function matroid_realization_space(n,nonbases)  
	matroid_variety = Variety(matroid_space_eqns(n, nonbases))
	NID = nid(matroid_variety)
	sorted_nonbases = sort([sort(m) for m in nonbases])
	WS = witness_sets(NID)

	CorrectComponents = []

	for i in keys(WS) #Scroll through each dimension

		number_components_dim_i = length(WS[i]) #number of components with dimension i 
		
		for j in 1:number_components_dim_i #Scroll through each of the components of that dimension

			solution = sample(WS[i][j]) #sample a point on that component
			matrix_solution = solution_to_matrix(solution) #turn that point into a matrix (projective)
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



#This takes a (complex) witness set and homotopies it to the linear subspace of just the REAL parts of the matrix A and vector b from which the witness set corresponded to the linear space L: Ax=b
function real_witness_set(W)
           newW = witness_set(W,LinearSubspace(real(W.L.extrinsic.A),real(W.L.extrinsic.b)))
           realWitness = HomotopyContinuation.real_solutions(newW.R)
           return realWitness
       end
       
function visualize_matroid(n, nonbases) 

	witness = decompose_matroid_space(n, nonbases) 
	pointVector = real_witness_set(witness) 	#real points to plot
	pointMatrix = reshape(pointVector[1], (2,n)) 
	x = pointMatrix[1,:]
	y = pointMatrix[2,:]
	visualization = scatter(x, y, legend = false)
		for i in 1:n
    			annotate!(x[i].+0.05, y[i].+0.05, text("$i", :left)) #labelling each point by number
		end
	
		for j in 1:length(nonbases)  #plotting lines the points lie on
			linelab = nonbases[j] 
			
			plot!([x[k] for k in nonbases[j]], [y[k] for k in nonbases[j]])
		end 
	return(visualization)
	
	#scores 
	point_pairs = collect(combinations(1:n, 2))
	points = []
	distance = []
	
		for k in 1:n
		push!(points, pointMatrix[:,k])
		end
	
	for j in 1:length(point_pairs) 		
		u = points[(point_pairs[j][1])]
		v = points[(point_pairs[j][2])]
		push!(distance, distance(u, v, EuclideanNorm()))
	end 						
	
end 
