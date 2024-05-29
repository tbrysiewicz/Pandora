export
	visualization_with_refinement
	

	

	

#TODO: Catch when the user tries to input an EP w/ more than two params
#      and suggest restricting


#Input: 
#		EP - an enumerative problem 
#		P  - a collection of (n) parameter values of EP
#Output:
#       a new enumerative problem which is EP restricted to the affine (n-1)-dimensional space
#       spanned by the parameters in P

function restrict_enumerative_problem(EP::EnumerativeProblem,P::Vector{Vector{Float64}})
	F = system(EP)
	xv = variables(F)
	xp = parameters(F)
	n = length(P)
	@var t[1:n-1]
	affine_span = P[n]+sum([t[i].*(P[i]-P[n]) for i in 1:n-1])
	NewEquations = [subs(f,xp=>affine_span) for f in expressions(F)]
	return(EnumerativeProblem(System(NewEquations,variables=xv,parameters=t)))
end
#=
EP = TwentySevenLines()
MyPlots = visualizationWithRefinement(restrict_enumerative_problem(EP,[randn(Float64,20) for i in 1:3]),[-1,1],[-1,1],100,5)
=#


#Input:
#		xlims - x parameter interval
#		ylims - y parameter interval
#		resolution - the maximum number of points outputted in the mesh	
#Output:
#		A collection of target parameters "on a grid", restricted by the inputted xlims and ylims with maximum number of outputted parameters limited by resoluton.

function create_mesh(xlims::Vector, ylims::Vector, resolution::Int)
	xRange = xlims[2]-xlims[1]
	yRange = ylims[2]-ylims[1]
	
	numberOfDivisions = floor(sqrt(resolution))-1
	
	xDivisionSize = xRange/numberOfDivisions
	yDivisionSize = yRange/numberOfDivisions
	
	targetParameters = [[a,b] for a in xlims[1]:xDivisionSize:xlims[2] for b in ylims[1]:yDivisionSize:ylims[2]]
	
	return targetParameters
end

#Input:
#		F - a system
#		data1 - a dataset produced by solving F using the "solve" function
#		numberOfSolutions - the "correct" number of solutions to a problem, i.e., the degree of the problem
#		certification - option to certify the solutions in data1
#Output:
#		A dictionary that takes each of the target parameters from data1 as keys and assigns the corresponding number of real solutions to its value.

function data_to_dictionary_converter(F::System, data1::Vector{Tuple{Result, Vector{Float64}}}, numberOfSolutions::Int, certification = false)
	D = Dict()
	for i in 1:length(data1)
		if certification == true
			C = certify(F, data1[i][1], data1[i][2])
				if length(findall(cert->is_certified(cert)==false, C.certificates))!=0 || length(solutions(data1[i][1]))!=numberOfSolutions || nsingular(data1[i][1])!=0
					D[data1[i][2]] = -2
				else
					D[data1[i][2]] = nreal(data1[i][1])
				end
		elseif certification == false
			if length(solutions(data1[i][1]))!=numberOfSolutions || nsingular(data1[i][1])!= 0
				D[data1[i][2]] = -2
			
			else
				D[data1[i][2]] = nreal(data1[i][1])	
			end
		end
		
	end
	return D
end

#Input:
#		parameterVector1 - a two-dimensional parameter
#		parameterVector2 - a two-dimensional parameter
#Output:
#		a collection of parameters on the "box" produced by the two inputted parameter pairs. If the pairs are colinear, the outputted parameter is simply their midpoint. Otherwise, seven parameters are returned
function box_refinement(parameterVector1::Vector, parameterVector2::Vector)
	boxParameters = []
		if parameterVector1[1] == parameterVector2[1] || parameterVector1[2]==parameterVector2[2]
			midpointVector = (0.5)*(parameterVector2 - parameterVector1) + parameterVector1
			push!(boxParameters,midpointVector)
		else 
			boxParameter1 = [parameterVector2[1], parameterVector1[2]]
			boxParameter2 = [parameterVector1[1], parameterVector2[2]]
			boxParameter3 = (0.5)*(boxParameter1-parameterVector1)+parameterVector1
			boxParameter4 = (0.5)*(parameterVector2-boxParameter2)+boxParameter2
			boxParameter5 = (0.5)*(boxParameter2 - parameterVector1) + parameterVector1
			boxParameter6 = (0.5)*(parameterVector2 - boxParameter1) + boxParameter1
			boxParameter7 = (0.5)*(boxParameter4 - boxParameter3) + boxParameter3
			
			push!(boxParameters, boxParameter1)
			push!(boxParameters, boxParameter2)
			push!(boxParameters, boxParameter3)
			push!(boxParameters, boxParameter4)
			push!(boxParameters, boxParameter5)
			push!(boxParameters, boxParameter6)
			push!(boxParameters, boxParameter7)
		end
	return boxParameters
end

#Input:
#		dictionary1 - a dictionary of the sort produced by data_to_dictionary_converter
#		xDistance - A numerical value that determines the x length of the rectangle used for refinement
#		yDistance - A numerical value that determines the y length of the rectangle used for refinement
#Output:
#		a collection of parameters. For each parameter in dictionary1, the function checks in the 2*xDistance x 2*yDistance rectangle centered at the parameter for other parameters with a differing number of real solutions. If such a parameter point exists, box_refinement is used to produce new target parameters. The function does this for every parameter in dictionary1 and then returns all of the generated parameters.


function refined_parameters(dictionary1::Dict, xDistance, yDistance)
	newParameters = []
	
	for i in keys(dictionary1)
		for j in keys(dictionary1)
			if  (abs(i[1]-j[1]))<=xDistance && (abs(i[2]-j[2]))<=yDistance && dictionary1[i]!=dictionary1[j]
				tempBoxParameters = box_refinement(i,j)
				
				for k in 1:length(tempBoxParameters)
					if (tempBoxParameters[k] in newParameters) == false
						push!(newParameters,tempBoxParameters[k])
					end
				end
			end
		end
	end
	
	return newParameters
end

#Input:
#		F - a system
#		dictionary1 - a dictionary of the sort produced by data_to_dictionary_converter
#		xDistance - A numerical value that determines the x length of the rectangle used for box_refinement
#		yDistance - A numerical value that determines the y length of the rectangle used for box_refinement
#		numberOfSolutions - the "correct" number of solutions to a problem, i.e., the degree of the problem
#		certification - option to certify solutions
#Output:
#		returns a dictionary with parameter keys with values corresponding to number of real solutions. The function takes dictionary1 and then uses refined_parameters to generate new parameters, solves for these parameters, and then returns a dictionary containing both these new points and the original points from dictionary1.

function refined_data(F::System, dictionary1::Dict, xDistance, yDistance, numberOfSolutions::Int, certification=False)
	newTargetParameters = refined_parameters(dictionary1, xDistance, yDistance)
	repeatParameters = [] 
	for i in 1:length(newTargetParameters)
		if (newTargetParameters[i] in keys(dictionary1)) == true
			push!(repeatParameters, newTargetParameters[i])
		end
	end
	setdiff!(newTargetParameters, repeatParameters)
	P = randn(ComplexF64, 2)
	S = solve(F; target_parameters = P)
	newData = solve(F, S; start_parameters = P, target_parameters = newTargetParameters)
	
	dictionary2 = data_to_dictionary_converter(F, newData, numberOfSolutions, certification)
	
	dictionary3 = merge(dictionary1,dictionary2)
	
	return dictionary3
end

#Input:
#		dictionary1 - a dictionary of the sort porduced by data_to_dictionary_converter
#Output:
#		a scatter plot of the parameters from dictionary1 with differing marker colors for different number of real solutions.

function parameter_dictionary_scatter(dictionary1)

	possibleNumberOfRealSolutions = []
	for i in keys(dictionary1)
		if (dictionary1[i] in possibleNumberOfRealSolutions) == false
			push!(possibleNumberOfRealSolutions, dictionary1[i])
		end
	end

	MyPlot = scatter()
	
	for i in 1:length(possibleNumberOfRealSolutions)
		number = possibleNumberOfRealSolutions[i]
		tempDataPoints = filter(x->dictionary1[x]==possibleNumberOfRealSolutions[i], keys(dictionary1))
		scatter!([x[1] for x in tempDataPoints], [x[2] for x in tempDataPoints], palette =:rainbow, markerstrokewidth = 0, markershape =:rect, markersize = 1.5, labels = "$(number) Real Solutions")
	end
	
	display(MyPlot)
	
	return MyPlot

end


#Input:
#		E - an EnumerativeProblem
#		xlims - the x paramater interval over which the system will be solved
#		ylims - the y parameter interval over which the system will be solved
#		initialResolution - the number of parameters to be solved for in the initial plot
#		depth - the number of times the plot is refined
#		certification - option to certify solutions
# Output:
#		Returns a collection of plots. The function solves for and plots the data for the system corresponding to E over the inputted xlims and ylims with a maximum number of data points equal to initialResolution. It then uses refined_data on this initial data to produce a second plot. This procedure repeats for the specified depth.

function visualization_with_refinement(E, xlims = [-2,2], ylims = [-2,2], initialResolution = 2000, depth = 3; certification = false)
	F = system(E)
	
	if length(parameters(F))!=2
		throw(ArgumentError("System does not consist of two parameters."))
	end
	
	numberOfSolutions = degree(E)
	
	mesh1 = create_mesh(xlims, ylims, initialResolution)
	xRange = xlims[2]-xlims[1]
	yRange = ylims[2]-ylims[1]
	
	numberOfDivisions = floor(sqrt(initialResolution))-1
	
	xDivisionSize = xRange/numberOfDivisions
	yDivisionSize = yRange/numberOfDivisions
	
	xDistance = xDivisionSize*(1.05)
	yDistance = yDivisionSize*(1.05)
	
	P = randn(ComplexF64, 2)
	S = solve(F; target_parameters = P)
	data1 = solve(F, S; start_parameters = P, target_parameters = mesh1)
	
	dictionary1 = data_to_dictionary_converter(F, data1, numberOfSolutions, certification)
	
	MyPlots = []
	push!(MyPlots,parameter_dictionary_scatter(dictionary1))
	
	dictionary2 = refined_data(F, dictionary1, xDistance, yDistance, numberOfSolutions, certification)
	
	push!(MyPlots,parameter_dictionary_scatter(dictionary2))
	
	for i in 2:depth
		xDistance2 = xDistance/2^(i-1)
		yDistance2 = yDistance/2^(i-1)
		dictionary3 = refined_data(F, dictionary2, xDistance2, yDistance2, numberOfSolutions, certification)
		push!(MyPlots,parameter_dictionary_scatter(dictionary3))
		
		dictionary2 = dictionary3
	end
	return(MyPlots)
end


