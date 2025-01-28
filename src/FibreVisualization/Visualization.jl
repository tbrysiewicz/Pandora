#I pulled the following three functions from main branch, only changed getters
function Gram_Schmidt(basis_vectors::Vector)
	orthonormal_basis = []
	for i in 1:length(basis_vectors)
		if i == 1
			new_basis_vector = basis_vectors[i]
			new_basis_vector = new_basis_vector/norm(new_basis_vector)
			push!(orthonormal_basis, new_basis_vector)
		else
			new_basis_vector = basis_vectors[i] - sum([(sum(basis_vectors[i].*orthonormal_basis[x])/sum(orthonormal_basis[x].*orthonormal_basis[x])).*orthonormal_basis[x] for x in 1:i-1])
			new_basis_vector = new_basis_vector/norm(new_basis_vector)
			push!(orthonormal_basis, new_basis_vector)
		end
	end
	return orthonormal_basis
end

function restrict_enumerative_problem(EP::EnumerativeProblem,P::Vector{Vector{Float64}})
    xv = variables(EP)
    xp = parameters(EP)
    n = length(P)
    @var t[1:n-1]
    basis_vectors = [(P[i]-P[n]) for i in 1:n-1]
    basis_vectors = Gram_Schmidt(basis_vectors)
    affine_span = P[n] + sum([t[i].*basis_vectors[i] for i in 1:n-1])
    new_expressions = [subs(f,xp=>affine_span) for f in expressions(system(EP))]
    return(EnumerativeProblem(System(new_expressions,variables=xv,parameters=t)))
end

function restrict_enumerative_problem_to_plane(EP::EnumerativeProblem)
	P = [randn(Float64,n_parameters(EP)) for i in 1:3]
	return(restrict_enumerative_problem(EP,P))
end

function initialize_subdivision(EP::EnumerativeProblem, xlims::Vector, ylims::Vector, fibre_function = x->n_real_solutions(x), resolution = 1000)

	xlength = xlims[2] - xlims[1]
	ylength = ylims[2] - ylims[1]

	if xlength == ylength
		num_x_divs = Int(floor(sqrt(resolution)))
		num_y_divs = Int(floor(sqrt(resolution)))
	else
		xlength_proportion = (xlength)/(xlength + ylength)
		num_x_divs = sqrt((xlength_proportion*resolution)/(1 - xlength_proportion))
		num_y_divs = resolution/num_x_divs
		num_x_divs = Int(floor(num_x_divs))
		num_y_divs = Int(floor(num_y_divs))
	end
			
	x_values = range(xlims[1], xlims[2], num_x_divs)
	y_values = range(ylims[1], ylims[2], num_y_divs)
	shift_amount = (x_values[2] - x_values[1])/2
	parameters_1 = [[i - isodd(findfirst(x->x==j, y_values))*shift_amount, j] for i in x_values for j in y_values]
	
	subdivision_solutions = solve(EP, parameters_1)
    function_cache = []
	
	for i in subdivision_solutions
		push!(function_cache, (parameters_1[findfirst(x->x==i, subdivision_solutions)], fibre_function(i)))
	end

	polygons = []
	
	parameter_array = reshape(parameters_1, length(x_values), length(y_values))
	for i in 1:length(x_values)-1
		for j in 1:length(y_values)-1
			tl = parameter_array[i,j]
			tr = parameter_array[i+1,j]
			bl = parameter_array[i, j+1]
			br = parameter_array[i+1, j+1]
			tl_index = findfirst(x->x[1] == tl, function_cache)
			tr_index = findfirst(x->x[1] == tr, function_cache)
			bl_index = findfirst(x->x[1] == bl, function_cache)
			br_index = findfirst(x->x[1] == br, function_cache)

			if isodd(i)
				push!(polygons, [tl_index, tr_index, bl_index])
				push!(polygons, [bl_index, tr_index, br_index])
			else
				push!(polygons, [tl_index, tr_index, br_index])
				push!(polygons, [bl_index, tl_index, br_index])
			end
		end
	end
	
	initial_graphmesh = GraphMesh(function_cache)
	initial_subdivision = Subdivision(initial_graphmesh,polygons)
	
	return initial_subdivision
end


function visualize(GM::GraphMesh)
    #Create plot, display plot, return plot
end


function visualize(EP::EnumerativeProblem; fibre_function = n_real_points,kwargs)
    #Create graph mesh
    #Refine graph mesh to user specs
    #Simplify unnecessary mesh points from the point of view of the plot 
    #Call visualize on graph mesh and return that output, all while cache-ing the graphmesh
end

