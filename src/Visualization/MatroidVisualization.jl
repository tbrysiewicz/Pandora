export 
	draw_matroid_representative

#Function: matroid_space_eqns  
#Input: M = matrix repersentative of matroid
#       nonbases = points we want collinear/ nonbases of matroid  
#Output: visualization = Point configuration of matroid in the plane 
function draw_matroid_representative(M::Matrix{Float64},nonbases::Vector{Vector{Int}})
	x = M[1,:] 
	y = M[2,:]
	visualization = scatter(x, y, markersize =4, markercolor=:black, legend = false) #scatter plot with columns of M as points
		for i in 1:size(M,2)
    			annotate!(x[i], y[i], text("$i", 4, :white, :center)) #labelling each point by number
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


#Function: plot_first_n_matroids- generates the first n simple matroids from the ploymake database that are realizable and do not have repeated nonbases, and plots them 
#input : n = interval of matroids to plot, N = number of steps taken by the optimizer for each matroid, l = layout of the combined plot
#output : combine_plot = a combined plot of all of the matroids, with layout l 

function  list_of_matroids(n:: Tuple{Int64, Int64}; n_trials = 10, n_samples=10)
	matroidsnb = [] #empty list of the nonbases of the matroids 
	matroid_repersentatives = []
	p = 1 #number of elements the matroid has 
	
	
	db = Polymake.Polydb.get_db();
	collection = db["Matroids.Small"]; #generating the database 
		
	while length(matroidsnb) < n[2]	#while the number of matroids added to list is less than n[2]
		query = Dict("RANK"=>3,"N_ELEMENTS"=>p,"SIMPLE"=>true) #query for generating all matroids with p elements from the database 
		results1 = Polymake.Polydb.find(collection, query) 
		oscar_matroids = [Matroid(pm) for pm in results1] #list of all matroids with p elements 
		
		for i in 2:length(oscar_matroids) #searching through all matroids in oscar_matroids and generating nonbases and output of best_realizable_matroid (skip the first matroid, as nb = [])
				   nb = nonbases(oscar_matroids[i])
				   matroid_matrix = best_realizable_matroid(oscar_matroids[i]; n_trials=n_trials, n_samples=n_samples) #outputs either a plot of oscar_matroid[i] or nothing if matroid is not realizable/ has no real solutions
			   
				   if !(nb in matroidsnb) && !isnothing(matroid_matrix)  #filtering out all matroids that have the same nonbases as another matroids in the list matroidsnb and filtering out matroids that are not realizable and have no real solutions
					   push!(matroidsnb, nb)

					   if length(matroidsnb) >= n[1]
					   		push!(matroid_repersentatives, (oscar_matroids[i], matroid_matrix))

					   		if length(matroidsnb) == n[2] #break if there is n matroids in the list 
						   		break 
					   		end
						end 
				   end 
		end
		p =p + 1 #updating to generate another query 
	end 
	return(matroid_repersentatives)
	end 

function draw_matroid_collection(S, l)

	plot_vector = Vector{Plots.Plot}() #empty vector of plots 
	for i in 1:length(S)
		p = draw_matroid_representative(S[i][2], S[i][1])
		push!(plot_vector, p)
	end
	combined_plot = plot(plot_vector..., layout = l, margins = -2*Plots.mm) #creating the combined plot of the n matroids 
	return(combined_plot)
end 