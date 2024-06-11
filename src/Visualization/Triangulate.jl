using Random

function visualize_with_triangles(EP;depth = 4, resolution=1000,xlims=[-2,2],ylims=[-2,2],fibre_function = x->HomotopyContinuation.nreal(x[1]))
    (V,T) = initial_triangular_mesh(EP;
                                    fibre_function = fibre_function, 
                                    xlims = xlims,
                                    ylims=ylims,
                                    resolution=resolution)
    for i in 2:depth
        (V,T) = refine_triangular_mesh(EP,V,T;
                                    fibre_function = fibre_function, 
                                    xlims = xlims,
                                    ylims=ylims,
                                    resolution=resolution)
    end

    myplot = draw_triangular_mesh(V,T)
end


function refine_triangular_mesh(EP,value_dict,Triangles; fibre_function = x->HomotopyContinuation.nreal(x[1]), xlims = [-2,2],ylims=[-2,2],resolution=1000)
    non_complete = filter(T->value_dict[T[1]]!=value_dict[T[2]] || value_dict[T[1]] != value_dict[T[3]],Triangles)
    if length(non_complete)>resolution 
        non_complete = unique(rand(non_complete,resolution))
    end
    new_triangles = setdiff(Triangles,non_complete)
    new_parameters=[]
    for T in non_complete
        barycenter = (T[1]+T[2]+T[3])/3
        push!(new_triangles,[T[1],T[2],barycenter])
        push!(new_triangles,[T[1],T[3],barycenter])
        push!(new_triangles,[T[2],T[3],barycenter])
        push!(new_parameters,barycenter)
    end
    S = solve_over_params(EP,new_parameters; checks = [])

    for s in S
        value_dict[s[2]]=fibre_function(s)
    end
    return((value_dict,new_triangles))
end

function initial_triangular_mesh(EP;fibre_function = x->HomotopyContinuation.nreal(x[1]), xlims = [-2,2],ylims=[-2,2],resolution=1000)
    nxy = Int(floor(sqrt(resolution+100)))
    xrange = range(xlims[1],xlims[2],nxy)
    yrange = range(ylims[1],ylims[2],nxy)
    Triangles = []
    P = [[i,j] for i in xrange for j in yrange]
    S = solve_over_params(EP,P; checks = [])
    value_dict=Dict{Any,Any}()
    for s in S
        value_dict[s[2]]=fibre_function(s)
    end
    for i in 1:length(xrange)-1
        for j in 1:length(yrange)-1
            bl = [xrange[i],yrange[j]]
            br = [xrange[i+1],yrange[j]]
            tl = [xrange[i],yrange[j+1]]
            tr = [xrange[i+1],yrange[j+1]]
            if rand([1])==1
                push!(Triangles,[tl,tr,bl])
                push!(Triangles,[bl,tr,br])
            else
                push!(Triangles,[tl,tr,br])
                push!(Triangles,[bl,tl,br])
            end
        end
    end
    return((value_dict,Triangles))
end

function draw_triangle!(T,v)
    Triangle = Shape([(t[1],t[2]) for t in T])
    c =cgrad(:thermal, rev = false)[v]
    plot!(Triangle, fillcolor = c, linecolor=c, linewidth=0)
end

function draw_triangular_mesh(value_dict,Triangles;xlims = [-2,2],ylims=[-2,2])
    M = max(values(value_dict)...)
    m = min(values(value_dict)...)
    myplot = plot(xlims=xlims,ylims=ylims,legend=false)
    for T in Triangles
        if value_dict[T[1]]==value_dict[T[2]] && value_dict[T[1]]==value_dict[T[3]]
            draw_triangle!(T,value_dict[T[1]]/M)
        end
    end
    return(myplot)
end


#Below is code for "triforce" triangulation visualization strategy


function refine_triangular_mesh2(EP::EnumerativeProblem, value_dictionary, triangles; fibre_function = x->HomotopyContinuation.nreal(x[1]), xlims = [-2,2],ylims=[-2,2],resolution=1000)
	whiteSpace = filter(x->value_dictionary[x[1]]!=value_dictionary[x[2]] || value_dictionary[x[1]]!=value_dictionary[x[3]], triangles)
	new_triangles = setdiff(triangles, whiteSpace)
	newParameters = []
	
	if (length(whiteSpace)*3)>resolution
		whiteSpace = unique(rand(whiteSpace,Int(floor(resolution/3))))
	end
	
		
	for i in whiteSpace
		midpoint1 = 0.5*(i[2]-i[1]) + i[1]
		midpoint2 = 0.5*(i[3]-i[2]) + i[2]
		midpoint3 = 0.5*(i[3]-i[1]) + i[1]
		
		push!(new_triangles, [midpoint1, midpoint3, i[1]])
		push!(new_triangles, [i[2], midpoint2, midpoint1])
		push!(new_triangles, [midpoint2, i[3], midpoint3])
		push!(new_triangles, [midpoint1, midpoint2, midpoint3])
		
		push!(newParameters, midpoint1)
		push!(newParameters, midpoint2)
		push!(newParameters, midpoint3)
	end
	
	S = solve_over_params(EP, newParameters, checks = [])
	
	for i in S
		value_dictionary[i[2]] = fibre_function(i)
	end
	
	return((value_dictionary, new_triangles))
end
		

function triforce_visualization(EP;depth = 4, resolution=1000,xlims=[-2,2],ylims=[-2,2],fibre_function = x->HomotopyContinuation.nreal(x[1]))
    (V,T) = initial_triangular_mesh(EP;
                                    fibre_function = fibre_function, 
                                    xlims = xlims,
                                    ylims=ylims,
                                    resolution=resolution)
    for i in 2:depth
        (V,T) = refine_triangular_mesh2(EP,V,T;
                                    fibre_function = fibre_function, 
                                    xlims = xlims,
                                    ylims=ylims,
                                    resolution=resolution)
    end

    myplot = draw_triangular_mesh(V,T)
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


