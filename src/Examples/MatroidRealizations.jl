


using HomotopyContinuation
using LinearAlgebra 
using Combinatorics 
using Plots

#Function: matroid_space_eqns  #Input: n = number of points, NonBases = points we want collinear/ nonbases of matroid  #Output: Eqns = determinant equations defining a matroid 
function matroid_space_eqns(n, nonbases) 
	@var x[1:2,1:n]
	matrix = vcat([1 for i in 1:n]', x) 
	eqns = [det(matrix[1:3,c]) for c in nonbases] #Computing equations 
	return eqns	
end 


#Function : solution_to_matrix #Input: soln = solution vector  #Output: m = matrix of a row of ones and solution
function solution_to_matrix(soln) 
	n = Int64(length(soln)/2) #Number of columns 
 	matrix_without_ones = reshape(soln, (2,n))  #Matrix of points without the row of ones 
 	matrix = vcat([1 for i in 1:n]' , matrix_without_ones) #Adding the row of 1s to make M
	return matrix
 end 


#Function : matrix_to_nonbases #Input: matrix = matirx   #Output: nonbases = Non-Bases of matrix
function matrix_to_nonbases(matrix) 
	nonbases= []
	n = size(matrix, 2) #number of columns
	r = collect(combinations(1:n, 3)) #collection of all possible bases 
	
	for i in 1:length(r)
		submatrix = matrix[:, r[i]] #All possible submatrices
		determinant_submatrix = det(submatrix) 
	
		if abs(determinant_submatrix) <0.000000001
			push!(nonbases, r[i]) 
		end	
	end
	return nonbases	
end 
 

#Function : decompose_matroid_space #Input: n= number of points, nonbases = points we want collinear/ nonbases of a matrix  #Output: matroidwitnessresults = Witness set of a matroid
function decompose_matroid_space(n,nonbases)  
	matroid_variety = Variety(System(matroid_space_eqns(n, nonbases)))
	witness_superset = nid(matroid_variety)
	sorted_nonbases = sort([sort(m) for m in nonbases])

	for i in keys(witness_superset.Witness_Sets)
		number_components_dim_i = length(witness_superset.Witness_Sets[i]) #number of components with dimension i 
		
		for j in 1:number_components_dim_i #Getting the nonbases of witness set components of dimension i 
			solution = sample(witness_superset.Witness_Sets[i][j]) #solution of witness set with dimension i 
			matrix_solution = solution_to_matrix(solution) #Turning this solution into matrix
			nonbases_component = matrix_to_nonbases(matrix_solution)  
			sorted_nonbases_component= sort([sort(m) for m in nonbases_component])
                 
			if sorted_nonbases_component == sorted_nonbases #Adding witness set component of dimension i to vector if non-bases are equal to input 
				return witness_superset.Witness_Sets[i][j]  
			end 	
		end 
		
	end 
end


function sample(W) #Takes witness set and samples random solution from it 
	return(rand(solutions(W)))
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
	
end 
