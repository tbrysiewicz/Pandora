function draw_matroid_representative(M::Matrix{Float64},nonbases::Vector{Vector{Int}})
	x = M[1,:]
	y = M[2,:]
	visualization = scatter(x, y, markersize =8, markercolor=:black, legend = false)
		for i in 1:size(M,2)
    			annotate!(x[i], y[i], text("$i", 8, :white, :center)) #labelling each point by number
		end
	
		for j in 1:length(nonbases)  #plotting lines the points lie on
			slope = (y[nonbases[j][2]] - y[nonbases[j][1]]) / (x[nonbases[j][2]] - x[nonbases[j][1]])
			b = y[nonbases[j][1]] - slope*x[nonbases[j][1]]
			f(x) = slope*x + b
			plot!(f)
		end 
	return(visualization)
	end

	function draw_matroid_representative(M::Matrix{Float64},Matr::Matroid)
		NB = nonbases(Matr)
		return(draw_matroid_representative(M,NB))
	end
