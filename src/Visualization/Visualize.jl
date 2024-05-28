export
	visualizationWithRefinement
	

	

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


#function Visualize(Something)
#	println("Probably delete this later. It is an example to learn how to use git");
#end

#Noah this is your job now.


#createMesh takes x value and y value limits as well as resolution and generates a grid of target parameter points.

function createMesh(xlims, ylims, resolution)
	xRange = xlims[2]-xlims[1]
	yRange = ylims[2]-ylims[1]
	
	numberOfDivisions = floor(sqrt(resolution))-1
	
	xDivisionSize = xRange/numberOfDivisions
	yDivisionSize = yRange/numberOfDivisions
	
	targetParameters = [[a,b] for a in xlims[1]:xDivisionSize:xlims[2] for b in ylims[1]:yDivisionSize:ylims[2]]
	
	return targetParameters
end

#dataDict takes a dataset from solving a polynomial system and produces a dictionary that takes parameter pairs as keys and number of real solutions as values. Option to certify solutions.

function dataDict(F, data1, numberOfSolutions, certification = true)
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

#boxRefinement takes two pairs of parameters and produces the additional parameters on the box produced by the two pairs. If the pairs are colinear then the additional parameter returned is their midpoint.

function boxRefinement(parameterVector1, parameterVector2)
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

#refinedParameters takes a dictionary of parameter pairs corresponding to a number of real solutions for those given parameters and goes through each key, checking if there are any other keys within a certain distance in terms of x or a certain distance in terms of y and has a different number of real solutions. If there are, it adds parameters produced by box refinement to a set of new parameters. These parameters are then returned and can be solved for.

function refinedParameters(dictionary1, xDistance, yDistance)
	newParameters = []
	
	for i in keys(dictionary1)
		for j in keys(dictionary1)
			if  (abs(i[1]-j[1]))<=xDistance && (abs(i[2]-j[2]))<=yDistance && dictionary1[i]!=dictionary1[j]
				tempBoxParameters = boxRefinement(i,j)
				
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

#refinedData uses refinedParameters to generate new parameters to solve for. It then solves for these parameters and returns a dictionary that contains all new points and number of real solutions as well as all of the initial points and their corresponding number of real solutions.

function refinedData(F, dictionary1, xDistance, yDistance, numberOfSolutions, certification)
	newTargetParameters = refinedParameters(dictionary1, xDistance, yDistance)
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
	
	dictionary2 = dataDict(F, newData, numberOfSolutions, certification)
	
	dictionary3 = merge(dictionary1,dictionary2)
	
	return dictionary3
end

#parameterDictionaryScatter takes as input a dictionary in which the keys are parameter points and values are number of real solutions and then plots a corresponding scatter plot.

function parameterDictionaryScatter(dictionary1)

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


#visualizationWithRefinement function will take a system and produce a plot of the parameter space with the given initial resolution. It will then refine the data and produce and save a plot at each refinement step.

function visualizationWithRefinement(E, xlims = [-5,5], ylims = [-5,5], initialResolution = 5000, refinementSteps = 3; certification = false)
	F = system(E)
	numberOfSolutions = degree(E)
	
	mesh1 = createMesh(xlims, ylims, initialResolution)
	xRange = xlims[2]-xlims[1]
	yRange = ylims[2]-ylims[1]
	
	numberOfDivisions = floor(sqrt(initialResolution))-1
	
	xDivisionSize = xRange/numberOfDivisions
	yDivisionSize = yRange/numberOfDivisions
	
	xDistance = xDivisionSize*(1.05)
	yDistance = yDivisionSize*(1.05)
	
	P = randn(ComplexF64, 2)
	S = solve(F; target_parameters = P)
	println(typeof(S))
	println(typeof(P))
	println(typeof(mesh1))
	data1 = solve(F, S; start_parameters = P, target_parameters = mesh1)
	
	dictionary1 = dataDict(F, data1, numberOfSolutions, certification)
	parameterDictionaryScatter(dictionary1)
	savefig("OriginalPlot.pdf")
	
	dictionary2 = refinedData(F, dictionary1, xDistance, yDistance, numberOfSolutions, certification)
	
	MyPlots = []
	push!(MyPlots,parameterDictionaryScatter(dictionary2))
	#savefig("Refinement1Plot.pdf")
	
	for i in 2:refinementSteps
		xDistance2 = xDistance/2^(i-1)
		yDistance2 = yDistance/2^(i-1)
		dictionary3 = refinedData(F, dictionary2, xDistance2, yDistance2, numberOfSolutions, certification)
		push!(MyPlots,parameterDictionaryScatter(dictionary3))
		savefig("Refinement$(i)Plot.pdf")
		
		dictionary2 = dictionary3
	end
	return(MyPlots)
end


