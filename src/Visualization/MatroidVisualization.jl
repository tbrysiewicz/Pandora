export 
	draw_matroid_representative
#LIST OF FUNCTIONS AND NOTES
#
#draw_matroid_repersentatives
#	- Check with Taylor if my documentation examples are ok? I can't really show the output since it is a plot.
#
#first_n_matroids
#first_n_simple_matroids
#first_n_simple_realizable_matroids
#	- I (Alex) added all of these functions on June 30th. They each do the same thing with different levels of filtering. I'm thinking I should just make one function, with options for filtering e.g. first_n_matroids(n::Int64 ; Simple = false, realizable = false)
#	- These functions output the first n matroids, but I could change it to output an interval of matroids, e.g. like matroids 5-10 from the database. 
#		
#draw_matroid_collection
#	-This was just to allow me to plot a layout of matroids together. Can probably be taken out, if we are going to plot matroids individually instead. 
#	
	@doc raw"""
	draw_matroid_representative(M::Union{Matrix{Int64}, Matrix{Float64}},nonbases::Vector{Vector{Int}})
	function draw_matroid_representative(M::Matrix{Float64},Matr::Matroid)
 
Returns a plot containing an affine drawing of a rank three matroid. $M$ is a matrix repersentative of the matroid described by $nonbases$ or $Matr$. The points of the plot correspond to the columns of the matrix $M$ and the nonbases are repersented by lines going through the corresponding points. 


# Examples
```jldoctest
julia> M = [1 0 1 4; 0 1 1 5; 0 0 0 6]
3×4 Matrix{Int64}:
1  0  1  4
0  1  1  5
0  0  0  6

julia> draw_matroid_representative(M, [[1,2,3]])
qt.qpa.plugin: Could not find the Qt platform plugin "wayland" in ""

julia> M = best_realizable_matroid(non_fano_matroid())
2×7 Matrix{Float64}:
0.0859759  -0.127336  -1.00552  1.40175    -0.761959   0.239214  0.480615
-0.518302   -0.204898   1.08535  0.0630736  -0.892962  -0.140661  0.454247

julia> draw_matroid_representative(M, non_fano_matroid())
qt.qpa.plugin: Could not find the Qt platform plugin "wayland" in ""

```
"""
function draw_matroid_representative(M::Union{Matrix{Int64},Matrix{Float64}},nonbases::Vector{Vector{Int}})
	x = M[1,:] 
	y = M[2,:]
	visualization = scatter(x, y, markersize =8, markercolor=:black, legend = false, grid = false) #scatter plot with columns of M as points
		for i in 1:size(M,2)
    			annotate!(x[i], y[i], text("$i", 8, :white, :center)) #labelling each point by number
		end
	
		for j in 1:length(nonbases)  #plotting lines the points lie on
			slope = (y[nonbases[j][2]] - y[nonbases[j][1]]) / (x[nonbases[j][2]] - x[nonbases[j][1]])
			b = y[nonbases[j][1]] - slope*x[nonbases[j][1]]
			f(x) = slope*x + b
			plot!(f, xaxis = false, yaxis = false, xlims = (minimum(x)-1 , maximum(x) + 1), ylims = (minimum(y)-1, maximum(y) + 1), linecolor=:blue, z_order=:back)
		end 
	return(visualization)
end 

function draw_matroid_representative(M::Matrix{Float64},Matr::Matroid)
	NB = nonbases(Matr)
	return(draw_matroid_representative(M,NB))
end


function first_n_matroids(n:: Int64)
	n_matroids = [] #list of oscar matroids 
	p = 1 #number of elements the matroid has 
	
	db = Polymake.Polydb.get_db();
	collection = db["Matroids.Small"]; #generating the database 
		
	while length(n_matroids) < n #while the number of matroids added to list is less than n[2]
		query = Dict("RANK"=>3,"N_ELEMENTS"=>p) #query for generating all matroids with p elements from the database 
		results1 = Polymake.Polydb.find(collection, query) 
		oscar_matroids = [Matroid(pm) for pm in results1] #list of all matroids with p elements 
		
		for i in 1:length(oscar_matroids) #searching through all matroids in oscar_matroids and generating nonbases
				   nb = nonbases(oscar_matroids[i])
				   elements_in_nb = unique(reduce(append!, nb, init = [])) #list of elements that appear in a nonbasis, without repeats 
				   
				   if length(elements_in_nb) == p #if all p elements appear in some nonbasis 
					   push!(n_matroids, oscar_matroids[i])
					   	if length(n_matroids) == n #break if there is n matroids in the list 
						   	break 
					   	end 
				   end 
		end
		p =p + 1 #updating to generate another query 
	end 
	return(n_matroids)
end 

 function first_n_simple_matroids(n:: Int64)
	n_matroids = [] #list of oscar matroids 
	p = 1 #number of elements the matroid has 
	
	db = Polymake.Polydb.get_db();
	collection = db["Matroids.Small"]; #generating the database 
		
	while length(n_matroids) < n #while the number of matroids added to list is less than n[2]
		query = Dict("RANK"=>3,"N_ELEMENTS"=>p, "SIMPLE" => true) #query for generating all matroids with p elements from the database 
		results1 = Polymake.Polydb.find(collection, query) 
		oscar_matroids = [Matroid(pm) for pm in results1] #list of all matroids with p elements 
		
		for i in 1:length(oscar_matroids) #searching through all matroids in oscar_matroids and generating nonbases
				   nb = nonbases(oscar_matroids[i])
				   elements_in_nb = unique(reduce(append!, nb, init = [])) #list of elements that appear in a nonbasis, without repeats 
				   
				   if length(elements_in_nb) == p #if all p elements appear in some nonbasis 
					   push!(n_matroids, oscar_matroids[i])
					   	if length(n_matroids) == n #break if there is n matroids in the list 
						   	break 
					   	end 
				   end 
		end
		p =p + 1 #updating to generate another query 
	end 
	return(n_matroids)
end 


function first_n_simple_realizable_matroids(n:: Int64) 
	n_matroids = [] #list of oscar matroids 
	p = 1 #number of elements the matroid has 
	
	db = Polymake.Polydb.get_db();
	collection = db["Matroids.Small"]; #generating the database 
		
	while length(n_matroids) < n #while the number of matroids added to list is less than n[2]
		query = Dict("RANK"=>3,"N_ELEMENTS"=>p, "SIMPLE" => true) #query for generating all matroids with p elements from the database 
		results1 = Polymake.Polydb.find(collection, query) 
		oscar_matroids = [Matroid(pm) for pm in results1] #list of all matroids with p elements 
		
		for i in 1:length(oscar_matroids) #searching through all matroids in oscar_matroids and generating nonbases
				   nb = nonbases(oscar_matroids[i])
				   elements_in_nb = unique(reduce(append!, nb, init = [])) #list of elements that appear in a nonbasis, without repeats 
				   
				   if length(elements_in_nb) == p #if all p elements appear in some nonbasis 
					   V = matroid_realization_space(oscar_matroids[i])
					   if !isnothing(V) #if oscar_matroid[i] is realizable 
					   	push!(n_matroids, oscar_matroids[i])
					   	if length(n_matroids) == n #break if there is n matroids in the list 
						   	break 
					   	end 
					   end 
				   end 
		end
		p =p + 1 #updating to generate another query 
	end 
	return(n_matroids)
end 


#Function: draw_matroid_collection 
#Input: S = list of tuples where each tuple contains a matrix repersentative and a matroid, l = layout of the combined plot
#Output: combined_plot = a combined plot of all of the matroids, with layout l 
function draw_matroid_collection(S :: Vector{Any}, l :: Tuple{Int64, Int64})

	plot_vector = Vector{Plots.Plot}() #empty vector of plots 
	for i in 1:length(S) #Scrolls through all matroids in S, draws them and then adds them to p
		p = draw_matroid_representative(S[i][2], S[i][1])
		push!(plot_vector, p)
	end
	combined_plot = plot(plot_vector..., layout = l, margins = -1*Plots.mm) #creating the combined plot of the n matroids 
	return(combined_plot)
end 