export
	visualization_with_refinement, restrict_enumerative_problem, moving_window_refinement
	

	

	

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

function create_mesh(xlims::Vector, ylims::Vector, resolution)
	if resolution < 8
		resolution = 8
	end
	
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

function parameter_dictionary_scatter(dictionary1::Dict)

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
		
		if number == -2
			scatter!([x[1] for x in tempDataPoints], [x[2] for x in tempDataPoints], palette =:rainbow, markerstrokewidth = 0, markershape =:rect, markersize = 1.5, labels = "Error")
		
		else
		scatter!([x[1] for x in tempDataPoints], [x[2] for x in tempDataPoints], palette =:rainbow, markerstrokewidth = 0, markershape =:rect, markersize = 1.5, labels = "$(number) Real Solutions")
		end
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

function visualization_with_refinement(E::EnumerativeProblem, xlims = [-2,2], ylims = [-2,2], initialResolution = 2000, depth::Int = 3; certification = false)
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


#Input:
#		xlims - interval for the x parameter over which the windows will be generated
#		ylims - interval for the y parameter over which the windows will be generated
#		windowxLength - the x parameter length of each window
#		windowyLength - the y parameter length of each window
#		coverage - the amount of overlap between adjacent windows, e.g. if coverage = 20, each window will overlap with 20% of the area of each adjacent window.
#Output:
#		A collection of vectors where each vector corresponds to a rectangular "window" defined by four coordinates.

function window_generator(xlims::Vector, ylims::Vector, windowxLength, windowyLength, coverage)
	upperLeftWindowCoords = [[A,B] for A in xlims[1]:(1-(0.01*coverage))*windowxLength:xlims[2] for B in ylims[2]:-(1-(0.01*coverage))*windowyLength:ylims[1]]
	upperLeftWindowCoords = filter(x->x[1]!=xlims[2] && x[2]!=ylims[1], upperLeftWindowCoords) #removing any coordinates on the lower or right boundary of the parameter space
	
	windowCoords = Vector[]
	
	for i in 1:length(upperLeftWindowCoords)
	
		tempWindowCoords = []
		push!(tempWindowCoords, upperLeftWindowCoords[i])
		
		tempLowerLeftWindowCoord = upperLeftWindowCoords[i] - [0,windowyLength]
		if tempLowerLeftWindowCoord[2]<ylims[1]
			tempLowerLeftWindowCoord[2] = ylims[1]
		end
		
		push!(tempWindowCoords, tempLowerLeftWindowCoord)
		
		tempUpperRightWindowCoord = upperLeftWindowCoords[i] + [windowxLength, 0]
		if tempUpperRightWindowCoord[1]>xlims[2]
			tempUpperRightWindowCoord[1] = xlims[2]
		end
		
		push!(tempWindowCoords, tempUpperRightWindowCoord)
		
		
		tempLowerRightWindowCoord = [tempUpperRightWindowCoord[1], tempLowerLeftWindowCoord[2]]
		push!(tempWindowCoords, tempLowerRightWindowCoord)
		
		push!(windowCoords, tempWindowCoords)
	end
	
	return windowCoords
end 

#Input:
#		dictionary1 - a dictionary of the sort produced by data_to_dictionary_converter
#		windows - a collection of windows as produced by window_generator
#		initialResolution - the resolution limit inputted for the initial plot
#Output:
#		returns a collection of parameters. The function takes the data from dictionary1 and "overlays" the windows onto the parameter space, checking in each window for parameter pairs corresponding to two distinct numbers of real solutions. If such a pair of points exists, the window is deemed "interesting" and a mesh of new parameters is generated over the area cut out by that window in the parameter space. The initial resolution is distributed equally amongst the interesting windows. These parameters are returned.

function refined_parameters2(dictionary1::Dict, windows::Vector{Vector}, initialResolution)
	
	interestingWindows = []
	newParameters = Set{Vector}()
	
	for i in windows
		solutionsInWindow = []
		for j in keys(dictionary1)
			if i[1][1]<=j[1] && j[1]<=i[3][1] && j[2]<=i[1][2] && i[2][2]<=j[2] && (dictionary1[j] in solutionsInWindow)==false
				push!(solutionsInWindow, dictionary1[j])
			end
			
			length(solutionsInWindow)>1 && break
		end
		
		if length(solutionsInWindow)>1
			push!(interestingWindows, i)
		end
	end
	
	windowResolution = initialResolution/length(interestingWindows)
	
	for i in interestingWindows
		xlims = [i[1][1], i[3][1]]
		ylims = [i[2][2], i[1][2]]
		
		mesh1 = create_mesh(xlims, ylims, windowResolution)
		
		for j in mesh1
			push!(newParameters, j)
		end
	end
	return newParameters
end

#Input
#		F - a System
#		dictionary1 - a dictionary of the sort produced by data_to_dictionary_converter
#		numberOfSolutions - The "correct" number of solutions, i.e., the degree of the problem
#		windows - a collection of windows as produced by window_generator
#		initialResoluton - the maximum resolution of the initial plot
#		certification - option to certify solutions
#Output
#		returns a dictionary of the sort produced by data_to_dictionary_converter. The function uses refined_parameters2 to create new target parameters. After solving for them it returns a dictionary that contains both these new points and their corresponding numbers of real solutions as well as the initial points from dictionary1.

function window_refinement_step(F::System, dictionary1::Dict, numberOfSolutions::Int, windows::Vector{Vector}, initialResolution, certification=false)
	targetParameters = refined_parameters2(dictionary1, windows, initialResolution)
	
	P = randn(ComplexF64, 2)
	S = solve(F; target_parameters = P)
	
	newSolutions = solve(F, S; start_parameters = P, target_parameters = targetParameters)

	dictionary2 = data_to_dictionary_converter(F, newSolutions, numberOfSolutions, certification)
	dictionary3 = merge(dictionary1, dictionary2)
	
	return dictionary3
end

#Input:
#		E - an EnumerativeProblem
#		xlims - interval for the x parameter
#		ylims - interval for the y parameter
#		initialResolution - the maximum resolution of the initial plot
#		depth - the number of refinement steps
#		certification - option to certify solutions
#		windowxLength - the x parameter length of each window
#		windowyLength - the y parameter length of each window
#		coverage - the amount of overlap between adjacent windows, e.g. if coverage = 20, each window will overlap with 20% of the area of each adjacent window.
#Output:
#		returns a collection of plots. The first is the initial plot and each subsequent plot is produced by plotting the dictionary given by window_refinement_step. In each stage of refinement, the dimensions of the windows are halved.

function moving_window_refinement(E::EnumerativeProblem, xlims::Vector, ylims::Vector, initialResolution, depth::Int, certification = false; windowxLength = ((xlims[2]-xlims[1])/10), windowyLength = ((ylims[2]-ylims[1])/10), coverage=20)
	F = system(E)
	
	if length(parameters(F))!=2
		throw(ArgumentError("System does not consist of two parameters."))
	end
	
	numberOfSolutions = degree(E)
	
	mesh1 = create_mesh(xlims, ylims, initialResolution)
	
	P = randn(ComplexF64,2)
	S = solve(F; target_parameters = P)
	data1 = solve(F,S; start_parameters = P, target_parameters = mesh1)
	
	dictionary1 = data_to_dictionary_converter(F, data1, numberOfSolutions, certification)
	
	myPlots = []
	
	push!(myPlots, parameter_dictionary_scatter(dictionary1))
	
	windows1 = window_generator(xlims, ylims, windowxLength, windowyLength, coverage)
	dictionary2 = window_refinement_step(F, dictionary1, numberOfSolutions, windows1, initialResolution, certification)
	
	push!(myPlots, parameter_dictionary_scatter(dictionary2))
	
	for i in 2:depth
		windowxLength1 = windowxLength/((2)^(i-1))
		windowyLength1 = windowyLength/((2)^(i-1))
		windows2 = window_generator(xlims, ylims, windowxLength1, windowyLength1, coverage)
		dictionary3 = window_refinement_step(F, dictionary2, numberOfSolutions, windows2, initialResolution, certification)
		push!(myPlots, parameter_dictionary_scatter(dictionary3))
		
		dictionary2 = dictionary3
	end
end





#rectangle_refinement is essentially the same as visualize_with_triangles except instead it uses rectangles. For each unresolved rectangle (rectangles with vertices that have different fibre function values) five new parameter points are generated subdividing the original rectangle into four new rectangles. Also has autofill function that generates an additional plot after plots for all depth levels are generated. After the last level of refinement, autofill subdivides each unresolved rectangle into four and colors each subrectangle with the color corresponding to the fibre function value for that vertex.

function rectangle_refinement(E::EnumerativeProblem; xlims = [-2,2], ylims = [-2,2], resolution = 1000, depth = 3, fibre_function = x->HomotopyContinuation.nreal(x[1]), fill = true)
	myPlots1 = []
	
	mesh1 = initial_rectangular_mesh(E, xlims, ylims, resolution, fibre_function)
	push!(myPlots1, plot_rectangular_mesh(mesh1[1], mesh1[2], xlims, ylims))
	
	for i in 2:depth
		mesh2 = refine_rectangular_mesh(E, mesh1[1], mesh1[2], resolution, fibre_function)
		push!(myPlots1, plot_rectangular_mesh(mesh2[1], mesh2[2], xlims, ylims))
		
		mesh1 = mesh2
	end
	
	if fill == true
		push!(myPlots1, auto_fill(mesh1[1], mesh1[2], xlims, ylims))
	end
	
	return myPlots1
end

function initial_rectangular_mesh(E::EnumerativeProblem, xlims = [-2,2], ylims = [-2,2], resolution = 1000, fibre_function = x->HomotopyContinuation.nreal(x[1]))
	divs = Int(floor(sqrt(resolution)))
	xvals = range(xlims[1], xlims[2], divs)
	yvals = range(ylims[1], ylims[2], divs)
	
	mesh1 = [[a,b] for a in xvals for b in yvals]
	
	data1 = solve_over_params(E, mesh1, checks = [])
	dictionary1 = Dict()
	numberOfSols = degree(E)
	
	
	for i in data1
		if length(solutions(i[1]))!=numberOfSols|| HomotopyContinuation.nsingular(i[1])!= 0
				dictionary1[i[2]] = -2
		else
		dictionary1[i[2]] = fibre_function(i)
		end
	end
	
	errorParams = filter(x->dictionary1[x] == -2, keys(dictionary1))
	
	if length(errorParams)!=0
		data2 = solve_over_params(E, errorParams, checks = [])
	
		for i in data2
			if length(solutions(i[1]))!=numberOfSols|| HomotopyContinuation.nsingular(i[1])!= 0
				dictionary1[i[2]] = -2
			else
				dictionary1[i[2]] = fibre_function(i)
			end
		end
	end
	
	rectangles = []
	
	for i in 1:length(xvals)-1
		for j in 1:length(yvals)-1
			bl = [xvals[i], yvals[j]]
			br = [xvals[i+1], yvals[j]]
			tl = [xvals[i], yvals[j+1]]
			tr = [xvals[i+1], yvals[j+1]]

			
			push!(rectangles, [bl, br, tr, tl])
		end
	end
	
	return (dictionary1, rectangles)
end

function plot_rectangle(rectangle, value)
	rectangle1 = Shape([(T[1], T[2]) for T in rectangle])
	c =cgrad(:thermal, rev = false)[value]
	plot!(rectangle1, fillcolor = c, linecolor=c, linewidth=0, labels = false)
end

function refine_rectangular_mesh(E::EnumerativeProblem, dictionary1, rectangles, resolution, fibre_function)

	numberOfSols = degree(E)
	
	whiteSpace = filter(i-> dictionary1[i[1]]!=dictionary1[i[2]] || dictionary1[i[1]]!=dictionary1[i[4]] || dictionary1[i[1]]!=dictionary1[i[3]], rectangles)
	
	if length(whiteSpace)*5 >resolution
		whiteSpace = unique(rand(whiteSpace,Int(floor(resolution/5))))
	end
	
	newRectangles = setdiff(rectangles, whiteSpace)
	
	newParams = []
	
	for i in whiteSpace
		center = [(((i[2][1]-i[1][1])/2)+i[1][1]), (((i[4][2]-i[1][2])/2)+i[1][2])]
		leftMidpoint = [i[1][1], center[2]]
		rightMidpoint = [i[2][1], center[2]]
		bottomMidpoint = [center[1], i[1][2]]
		topMidpoint = [center[1], i[3][2]]
		
		push!(newParams, center)
		push!(newParams, leftMidpoint)
		push!(newParams, rightMidpoint)
		push!(newParams, bottomMidpoint)
		push!(newParams, topMidpoint)
		
		push!(newRectangles, [i[1], bottomMidpoint, center, leftMidpoint])
		push!(newRectangles, [bottomMidpoint, i[2], rightMidpoint, center])
		push!(newRectangles, [center, rightMidpoint, i[3], topMidpoint])
		push!(newRectangles, [leftMidpoint, center, topMidpoint, i[4]]) 
	end
	
	data1 = solve_over_params(E, newParams, checks = [])
	dictionary2 = Dict()
	
	for i in data1
		if length(solutions(i[1]))!=numberOfSols|| nsingular(i[1])!= 0
				dictionary2[i[2]] = -2
		else
		dictionary2[i[2]] = fibre_function(i)
		end
	end
	
	errorParams = filter(x->dictionary2[x] == -2, keys(dictionary2))
	
	if length(errorParams)!=0
		data2 = solve_over_params(E, errorParams, checks = [])
	
		for i in data2
			if length(solutions(i[1]))!=numberOfSols|| nsingular(i[1])!= 0
				dictionary2[i[2]] = -2
			else
			dictionary2[i[2]] = fibre_function(i)
			end
		end
	end
	
	dictionary3 = merge(dictionary1, dictionary2)
	
	return(dictionary3, newRectangles)
end

function plot_rectangular_mesh(dictionary1, rectangles, xlims, ylims)
	myPlot = plot(xlims = xlims, ylims = ylims)
	M = max(values(dictionary1)...)
	
	for i in rectangles
		if dictionary1[i[1]]==dictionary1[i[2]]==dictionary1[i[3]]==dictionary1[i[4]]
			plot_rectangle(i, dictionary1[i[1]]/M)
		end
	end
	
	return myPlot
end


function auto_fill(dictionary1, rectangles, xlims, ylims)
	whiteSpace = filter(i-> dictionary1[i[1]]!=dictionary1[i[2]] || dictionary1[i[1]]!=dictionary1[i[4]] || dictionary1[i[1]]!=dictionary1[i[3]], rectangles)
	
	setdiff!(rectangles, whiteSpace)
	dictionary2 = Dict()
	
	for i in rectangles
		dictionary2[i] = dictionary1[i[1]]
	end
	
	for i in whiteSpace
		center = [(((i[2][1]-i[1][1])/2)+i[1][1]), (((i[4][2]-i[1][2])/2)+i[1][2])]
		leftMidpoint = [i[1][1], center[2]]
		rightMidpoint = [i[2][1], center[2]]
		bottomMidpoint = [center[1], i[1][2]]
		topMidpoint = [center[1], i[3][2]]
		
		bl = [i[1], bottomMidpoint, center, leftMidpoint]
		br = [bottomMidpoint, i[2], rightMidpoint, center]
		tr = [center, rightMidpoint, i[3], topMidpoint]
		tl = [leftMidpoint, center, topMidpoint, i[4]]
		
		dictionary2[bl] = dictionary1[i[1]]
		dictionary2[br] = dictionary1[i[2]]
		dictionary2[tr] = dictionary1[i[3]]
		dictionary2[tl] = dictionary1[i[4]]
		
	end
	
	myPlot1 = plot(xlims = xlims, ylims = ylims)
	M = max(values(dictionary2)...)
	
	for i in keys(dictionary2)
		plot_rectangle(i, dictionary2[i]/M)
	end
	
	println("Number of unresolved rectangles autofilled: ", length(whiteSpace))
	
	return myPlot1
end





