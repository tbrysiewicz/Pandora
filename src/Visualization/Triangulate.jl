#using Random

using DelaunayTriangulation

export
	visualize_discriminant

#General visualize_discriminant function that allows you to use different strategies: "Delaunay", "Triforce", "Barycenter"
function visualize_discriminant(EP::EnumerativeProblem, strategy; fibre_function = x->HomotopyContinuation.nreal(x[1]), xlims = [-2,2], ylims = [-2,2], resolution = 1000, depth = 4, automatic = false, total_resolution = 8*resolution, scatter = false, label = true, edges = false, continuous = false, plot_proportion = 0.7)
	if strategy == "Delaunay"
		return Delaunay_visualization(EP, xlims = xlims, ylims = ylims, fibre_function = fibre_function, depth = depth, resolution = resolution, total_resolution = total_resolution, automatic = automatic, scatter = scatter, edges = edges, label = label, continuous = continuous, plot_proportion = plot_proportion)
	elseif strategy == "Triforce"
		return triforce_visualization(EP, xlims = xlims, ylims = ylims, fibre_function = fibre_function, depth = depth, resolution = resolution, total_resolution = total_resolution, automatic = automatic, scatter = scatter, edges = edges, label = label, continuous = continuous, plot_proportion = plot_proportion)
	elseif strategy == "Barycenter"
		return Barycenter_visualization(EP, xlims = xlims, ylims = ylims, fibre_function = fibre_function, depth = depth, resolution = resolution, totalResolution = total_resolution, automatic = automatic, scatter = scatter, label = label, edges = edges, continuous = continuous, plot_proportion = plot_proportion)
	elseif strategy == "Delaunay2"
		return Delaunay_visualization_with_Ruppert_refinement(EP, xlims = xlims, ylims = ylims, fibre_function = fibre_function, depth = depth, resolution = resolution, total_resolution = total_resolution, automatic = automatic, scatter = scatter, label = label, edges = edges, continuous = continuous, plot_proportion = plot_proportion)
	else
		println("Invalid strategy inputted. Valid strategies include:")
		println("Delaunay")
		println("Delaunay2")
		println("Barycenter")
		println("Triforce")
	end
end	

#=function: Barycenter_visualization
Input:
		EP - EnumerativeProblem
		depth - number of refinement steps, used when automatic == false
		resolution - number of mesh vertices; number of parameters to be solved for
		totalResolution - total resolution used when automatic == true
		xlims - range of "x" parameter for mesh
		ylims - range of "y" parameter for mesh 
		fibre_function - function on solutions to EP
		automatic - option to run "automatic" refinement: refinement begins with a totalResolution and runs through each level of refinement using as much resolution as totalResolution permits.
		Refinement ends when totalResolution runs out. When automatic == false, the mesh is refined n = depth times with the same upper limit on number of refinement points = resolution.
		scatter - option to plot triangle vertices
		label - option to include a legend in plot
		edges - option to plot triangle edges 
		continuous - indicate whether fibre_function is discrete or continuous
		plot_proportion - in the case that fibre_function is continuous, plot_proportion is used to calculate the threshold for both refinement and plotting. For continuous fibre functions, the
		range of fibre_function values across a triangle is calculated (range across the values of the three vertices) and this range is used to determine if the triangle will be refined. For example,
		if plot_proportion is 0.7, the 70th percentile of triangle value ranges is used as the threshold for refinement; any triangle with a value range greater than the 70th percentile will be refined. 
		When refinement is complete, any triangle with a value range less than the 70th percentile will be plotted.
Output:
		myplot - A plot of the mesh following all stages of refinement
=#
function Barycenter_visualization(EP;depth = 4, resolution=1000, totalResolution = 8*resolution, xlims=[-2,2],ylims=[-2,2],fibre_function = x->HomotopyContinuation.nreal(x[1]), automatic = true, scatter = false, label = true, edges = false, continuous = false, plot_proportion = 0.7)
    V, T, L = initial_triangular_mesh(EP; fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = resolution)
	if automatic == false
		if continuous == true
			refinement_standard = standard_for_continuous_function_refinement(V, T, plot_proportion)
			for i in 2:depth
				V, T, A, L = refine_triangular_mesh(EP, V, T; fibre_function = fibre_function, xlims = xlims, ylims=ylims, resolution=resolution, continuous = true, standard_for_refinement = refinement_standard)
			end
		else
			for i in 2:depth
				V, T, A, L = refine_triangular_mesh(EP, V, T; fibre_function = fibre_function, xlims = xlims, ylims=ylims, resolution=resolution, continuous = false)
			end
		end
	else
		totalResolution -= L
		if continuous == true
			refinement_standard = standard_for_continuous_function_refinement(V, T, plot_proportion)
			while totalResolution > 0
				V, T, A, L = refine_triangular_mesh(EP, V, T; fibre_function = fibre_function, xlims = xlims, ylims=ylims, resolution = totalResolution, continuous = true, standard_for_refinement = refinement_standard)
				if A == true
					break
				end
				totalResolution -= L
			end
		else
			while totalResolution > 0
				V, T, A, L = refine_triangular_mesh(EP, V, T; fibre_function = fibre_function, xlims = xlims, ylims=ylims, resolution = totalResolution, continuous = false)
				if A == true
					break
				end
				totalResolution -= L
			end
		end
	end
	if continuous == true
		myplot = draw_triangular_mesh(V,T, scatter1 = scatter, label = label, edges = edges, continuous = true, plotting_standard = refinement_standard)
	else
    	myplot = draw_triangular_mesh(V,T, scatter1 = scatter, label = label, edges = edges)
	end
	return myplot
end

#=function: initial_triangular_mesh
input:
		EP - EnumerativeProblem
		fibre_function - function on solutions to EP
		xlims - range of "x" parameter for mesh
		ylims - range of "y" parameter for mesh 
		resolution - number of mesh vertices; number of parameters to be solved for
output:
		value_dict - Dictionary that takes parameters as keys and fibre_function output for values
		Triangles - Vector containing all triangles (defined by their vertices) that the mesh is composed of
=#
function initial_triangular_mesh(EP::EnumerativeProblem;fibre_function = x->HomotopyContinuation.nreal(x[1]), xlims = [-2,2],ylims=[-2,2],resolution=1000)
    nxy = Int(floor(sqrt(resolution))) #mesh points are distributed equally along x and y axis regardless of size of xlims and ylims producing "square" mesh
    xrange = range(xlims[1],xlims[2],nxy)
    yrange = range(ylims[1],ylims[2],nxy)
    delta = (xrange[2]-xrange[1])/2
    Triangles = []
    P = [[i-isodd(findfirst(x->x==j,yrange))*delta,j] for i in xrange for j in yrange]
    S = solve_over_params(EP,P; checks = [])
    value_dict=Dict{Any,Any}()
    for s in S
        value_dict[s[2]]=fibre_function(s)
        if length(solutions(s[1]))!=degree(EP) || nsingular(s[1])!=0
        	value_dict[s[2]] = false #parameters are given false value if they produce an "error"
        end
    end
    M = reshape(P,length(xrange),length(yrange))
    for i in 1:length(xrange)-1
        for j in 1:length(yrange)-1
            tl = M[i,j]#[xrange[i],yrange[j]+isodd(j)*(1/(2*nxy))]
            tr = M[i+1,j]#[xrange[i+1],yrange[j]+isodd(j)*(1/(2*nxy))]
            bl = M[i,j+1]#[xrange[i],yrange[j+1]+isodd(j)*(1/(2*nxy))]
            br = M[i+1,j+1]#[xrange[i+1],yrange[j+1]+isodd(j)*(1/(2*nxy))]
            if isodd(i)
                push!(Triangles,[tl,tr,bl])
                push!(Triangles,[bl,tr,br])
            else
                push!(Triangles,[tl,tr,br])
                push!(Triangles,[bl,tl,br])
            end
        end
    end
    return((value_dict, Triangles, length(P)))
end

#= function: refine_triangular_mesh
input:
		EP - EnumerativeProblem
		value_dict - Dictionary that takes parameters as keys and fibre_function output for values
		Triangles - Vector containing all triangles (defined by their vertices) that the initial mesh is composed of
		fibre_function - function on solutions to EP
		xlims - range of "x" parameter for mesh
		ylims - range of "y" parameter for mesh 
		resolution - number of parameters to be solved for during refinement; upper limit on number of triangles to be refined
		continuous - indicate whether fibre_function is discrete or continuous
		standard_for_refinement - in the case that fibre_function is continuous, used as the threshold for triangle 
		refinement; this value is produced by standard_for_continuous_function_refinement
output:
		value_dict - contains parameters from initial mesh as well as new parameters generated by refinement
		newTriangles - vector containing all triangles in refined mesh including triangles newly generated by refinement
		automaticTermination - used to terminate refinement loop in the case of automatic refinement
		length(S) - used to track parameters solved for in the case of automatic refinement
=#
function refine_triangular_mesh(EP::EnumerativeProblem,value_dict::Dict,Triangles; fibre_function = x->HomotopyContinuation.nreal(x[1]), xlims = [-2,2],ylims=[-2,2],resolution=1000, continuous = false, standard_for_refinement = 0.0)
	automaticTermination = false
	if continuous == true
		triangle_value_ranges = Dict()
        for T in Triangles
            min_value = min([value_dict[x] for x in T]...)
            max_value = max([value_dict[x] for x in T]...)
            triangle_value_ranges[T] = max_value - min_value
        end
        whiteSpace = filter(x->triangle_value_ranges[x]>standard_for_refinement, Triangles)
    else
        whiteSpace = filter(x->value_dict[x[1]]!=value_dict[x[2]]||value_dict[x[1]]!=value_dict[x[3]], triangles)
    end
	#whiteSpace contains all "unresolved" triangles: triangles with vertices that have different fibre_function values
    TrianglesToRefine=[]
	if length(whiteSpace)>resolution
		while length(TrianglesToRefine)<resolution #if whiteSpace contains more triangles than the refinement resolution, n = resolution triangles are randomly selected and later refined
			v = rand(whiteSpace,Int(floor(resolution/10)))
			TrianglesToRefine = unique(vcat(TrianglesToRefine,v))
		end
		TrianglesToRefine = TrianglesToRefine[1:Int(resolution)]
		automaticTermination = true 
	else
		TrianglesToRefine = whiteSpace
	end
    newTriangles = setdiff(Triangles, TrianglesToRefine) #newTriangles stores the trianlges that are not going to be refined
    newParameters=[]
    for T in TrianglesToRefine
        barycenter = (T[1]+T[2]+T[3])/3
        push!(newTriangles,[T[1],T[2],barycenter])
        push!(newTriangles,[T[1],T[3],barycenter])
        push!(newTriangles,[T[2],T[3],barycenter])
        push!(newParameters,barycenter)
    end
    S = solve_over_params(EP,newParameters; checks = [])
    for s in S
        value_dict[s[2]]=fibre_function(s)
        if length(solutions(s[1]))!=degree(EP) || nsingular(s[1])!=0
        	value_dict[s[2]] = false
        end
    end
    return(value_dict, newTriangles, automaticTermination, length(S))
end

#= function: draw_triangle
input:
		T - a vector containing a triangle (its vertices)
		v - numerical value to produce triangle's color
		edges - option to plot triangle with visible edges
		label - option to plot triangle with a corresponding label
		labelText - text to be written for label if label == true
output:
		no return. draw_triangle plots the inputted triangle
=#
function draw_triangle!(T, v; edges = false, label = false, labelText = "hello")
    Triangle = Shape([(t[1],t[2]) for t in T])
    c =cgrad(:thermal, rev = false)[v]
	if edges == false
		if label == true
			plot!(Triangle, fillcolor = c, linecolor=c, linewidth=0, labels = labelText)
		else
			plot!(Triangle, fillcolor = c, linecolor=c, linewidth=0, labels = false)
		end
	else
		if label == true
			plot!(Triangle, fillcolor = c, linecolor=:black, linewidth=0.2, labels = labelText)
		else
			plot!(Triangle, fillcolor = c, linecolor=:black, linewidth=0.2, labels = false)
		end
	end
end

#= function: draw_triangular_mesh
input:
		value_dict - Dictionary that takes parameters as keys and fibre_function output for values
		Triangles - Vector containing all triangles (defined by their vertices) that the final mesh is composed of
		xlims - range of "x" parameter for mesh
		ylims - range of "y" parameter for mesh
		scatter1 - option to plot triangle vertices
		label - option to include a legend in plot
		edges - option to plot triangle edges 
		continuous - indicate whether fibre_function is discrete or continuous
		plotting_standard - in the case that fibre_function is continuous plotting_standard is the threshold for triangle value ranges
		used to determine whether a triangle is to be plotted. This value is produced by standard_for_continuous_function_refinement
output:
		returns plot of parameter space
=#
function draw_triangular_mesh(value_dict::Dict,Triangles::Vector;xlims = [-2,2],ylims=[-2,2], continuous = false, scatter1 = true, label = true, edges = false, plotting_standard = 0.0)
    if continuous == true
        myplot = plot(xlims=xlims, ylims=ylims, colorbar = true)
        triangle_value_ranges = Dict()
        vertices = []
        for T in Triangles
            max_value = max([value_dict[x] for x in T]...)
            min_value = min([value_dict[x] for x in T]...)
            triangle_value_ranges[T] = max_value - min_value
        end
        trianglesToPlot = filter(x->triangle_value_ranges[x]<=plotting_standard, Triangles)
        triangle_plotting_values = Dict()
        for T in trianglesToPlot
            triangle_plotting_values[T] = Statistics.mean([value_dict[x] for x in T])
            for x in T
                push!(vertices, x)
            end
        end
        plotting_max = max([value_dict[x] for x in vertices]...)
        plotting_min = min([value_dict[x] for x in vertices]...)
        if scatter1 == true
            scatter!([a[1] for a in vertices], [a[2] for a in vertices], zcolor = [value_dict[a] for a in vertices], color =:thermal, markershape =:rect, markerstrokewidth = 0, labels = false, markersize = 0.2)
        else
            scatter!([a[1] for a in vertices], [a[2] for a in vertices], zcolor = [value_dict[a] for a in vertices], color =:thermal, markershape =:rect, markerstrokewidth = 0, labels = false, markersize = 0)
        end

        for T in trianglesToPlot
            shape1 = Shape([(t[1],t[2]) for t in T])
            color_value = (triangle_plotting_values[T]-plotting_min)/(plotting_max-plotting_min)
            if edges == true
                plot!(shape1, fillcolor = cgrad(:thermal, rev = false)[color_value], linecolor =:black, linewidth = 0.2, labels = false)
            else
                plot!(shape1, fillcolor = cgrad(:thermal, rev = false)[color_value], linecolor = cgrad(:thermal, rev = false)[color_value], linewidth = 0, labels = false)
            end
        end
    else
        M = max(values(value_dict)...)
        myplot = plot(xlims=xlims,ylims=ylims,legend=true)
        trianglesToPlot = filter(x->value_dict[x[1]]==value_dict[x[2]]&&value_dict[x[1]]==value_dict[x[3]], Triangles)
        if label == true
            valuesPlotted = []
            for i in trianglesToPlot
                if (value_dict[i[1]] in valuesPlotted) == false
                    if edges == false
                        draw_triangle!(i, value_dict[i[1]]/M, label = true, labelText = "$(value_dict[i[1]]) real solutions", edges = false)
                    else
                        draw_triangle!(i, value_dict[i[1]]/M, label = true, labelText = "$(value_dict[i[1]]) real solutions", edges = true)
                    end
                    push!(valuesPlotted, value_dict[i[1]])
                else
                    if edges == false
                        draw_triangle!(i, value_dict[i[1]]/M, label = false, edges = false)
                    else
                        draw_triangle!(i, value_dict[i[1]]/M, edges = true)
                    end
                end
            end
        else
            for i in trianglesToPlot
                if edges == false
                    draw_triangle!(i, value_dict[i[1]]/M, edges = false, label = false)
                else
                    draw_triangle!(i, value_dict[i[1]]/M, edges = true, label = false)
                end
            end
        end
        if scatter1 == true
            for i in unique(values(value_dict))
                temp_parameters = filter(x->value_dict[x]==i, keys(value_dict))
                if i == -2
                    scatter!([A[1] for A in temp_parameters], [A[2] for A in temp_parameters], markercolor =:red, markersize = 0.8, markershape =:rect, markerstrokewidth = 0.2, labels = false)
                else
                    c =cgrad(:thermal, rev = false)[i/M]
                    scatter!([A[1] for A in temp_parameters], [A[2] for A in temp_parameters], markercolor = c, markersize = 0.8, markershape =:rect, markerstrokewidth = 0.2, labels = false)
                end
            end
        end
    end
        return(myplot)
end

#=function: triforce_refinement 
input:
		EP - EnumerativeProblem
		value_dict - Dictionary that takes parameters as keys and fibre_function output for values
		triangles - Vector containing all triangles (defined by their vertices) that the initial mesh is composed of
		fibre_function - function on solutions to EP
		xlims - range of "x" parameter for mesh
		ylims - range of "y" parameter for mesh 
		resolution - number of parameters to be solved for during refinement; upper limit on number of triangles to be refined = floor(resolution/3)
		continuous - indicate whether fibre_function is discrete or continuous
		standard_for_refinement - in the case that fibre_function is continuous, used as the threshold for triangle 
		refinement; this value is produced by standard_for_continuous_function_refinement
output:
		value_dict - contains parameters from initial mesh as well as new parameters generated by refinement
		newTriangles - vector containing all triangles in refined mesh including triangles newly generated by refinement
		automaticTermination - used to terminate refinement loop in the case of automatic refinement
		length(S) - used to track parameters solved for in the case of automatic refinement
=#
function triforce_refinement(EP::EnumerativeProblem, value_dict::Dict, triangles::Vector; fibre_function = x->HomotopyContinuation.nreal(x[1]), xlims = [-2,2],ylims=[-2,2],resolution=1000, continuous = false, standard_for_refinement = 0.0)
	automaticTermination = false
    if continuous == true
        triangle_value_ranges = Dict()
        for T in triangles
            min_value = min([value_dict[x] for x in T]...)
            max_value = max([value_dict[x] for x in T]...)
            triangle_value_ranges[T] = max_value - min_value
        end
        whiteSpace = filter(x->triangle_value_ranges[x]>standard_for_refinement, triangles)
    else
        whiteSpace = filter(x->value_dict[x[1]]!=value_dict[x[2]]||value_dict[x[1]]!=value_dict[x[3]], triangles)
    end
    newParameters = []
	TrianglesToRefine=[]
	if (length(whiteSpace)*3)>resolution
		while length(TrianglesToRefine)*3<resolution
			v = rand(whiteSpace,Int(floor(resolution/10)))
			TrianglesToRefine = unique(vcat(TrianglesToRefine,v))
		end
		TrianglesToRefine = TrianglesToRefine[1:Int(floor(resolution/3))]
		automaticTermination = true
	else
		TrianglesToRefine = whiteSpace
	end
	newTriangles = setdiff(triangles, TrianglesToRefine) #Triangles not to refine. 
	for i in TrianglesToRefine
		midpoint1 = 0.5*(i[2]-i[1]) + i[1]
		midpoint2 = 0.5*(i[3]-i[2]) + i[2]
		midpoint3 = 0.5*(i[3]-i[1]) + i[1]
		push!(newTriangles, [midpoint1, midpoint3, i[1]])
		push!(newTriangles, [i[2], midpoint2, midpoint1])
		push!(newTriangles, [midpoint2, i[3], midpoint3])
		push!(newTriangles, [midpoint1, midpoint2, midpoint3])
		push!(newParameters, midpoint1)
		push!(newParameters, midpoint2)
		push!(newParameters, midpoint3)
	end
	S = solve_over_params(EP, newParameters, checks = [])
	for i in S
		if length(solutions(i[1]))!=degree(EP) || nsingular(i[1])!=0
        	value_dict[i[2]] = false
		else
			value_dict[i[2]] = fibre_function(i)
        end
	end
    return(value_dict, newTriangles, automaticTermination, length(S))
end

#=function: triforce_visualization
input:
		EP - EnumerativeProblem
		depth - number of refinement steps, used when automatic == false
		resolution - number of mesh vertices; number of parameters to be solved for
		total_resolution - total resolution used when automatic == true
		xlims - range of "x" parameter for mesh
		ylims - range of "y" parameter for mesh 
		fibre_function - function on solutions to EP
		automatic - option to run "automatic" refinement: refinement begins with a totalResolution and runs through each level of refinement using as much resolution as totalResolution permits.
		Refinement ends when totalResolution runs out. When automatic == false, the mesh is refined n = depth times with the same upper limit on number of refinement points = resolution.
		scatter - option to plot triangle vertices
		label - option to include a legend in plot
		edges - option to plot triangle edges 
		continuous - indicate whether fibre_function is discrete or continuous
		plot_proportion - in the case that fibre_function is continuous, plot_proportion is used to calculate the threshold for both refinement and plotting. For continuous fibre functions, the
		range of fibre_function values across a triangle is calculated (range across the values of the three vertices) and this range is used to determine if the triangle will be refined. For example,
		if plot_proportion is 0.7, the 70th percentile of triangle value ranges is used as the threshold for refinement; any triangle with a value range greater than the 70th percentile will be refined. 
		When refinement is complete, any triangle with a value range less than the 70th percentile will be plotted.

output:
		myplot - A plot of the mesh following all stages of refinement
=#
function triforce_visualization(EP;depth = 4, resolution=1000, total_resolution = 8*resolution, xlims=[-2,2],ylims=[-2,2],fibre_function = x->HomotopyContinuation.nreal(x[1]), continuous = false, scatter = false, automatic = true, edges = false, label = true, plot_proportion = 0.7)
	V, T, L = initial_triangular_mesh(EP, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = resolution)
	if automatic == true
		total_resolution -= L
        if continuous == true
            refinement_standard = standard_for_continuous_function_refinement(V,T, plot_proportion)
            while total_resolution > 0
                V, T, A, L= triforce_refinement(EP, V, T, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = total_resolution, continuous = true, standard_for_refinement = refinement_standard)
                if A == true
                    break
                end
                total_resolution -= L
		    end
        else
            while total_resolution > 0
                V, T, A, L = triforce_refinement(EP, V, T, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = total_resolution, continuous = false)
                if A == true
                    break
                end
                total_resolution -= L
            end
        end
	else
        if continuous == true
            refinement_standard = standard_for_continuous_function_refinement(V,T, plot_proportion)
            for i in 2:depth
                V, T, A, L = triforce_refinement(EP, V, T, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = resolution, continuous = true, standard_for_refinement = refinement_standard)
            end
        else
            for i in 2:depth
                V, T, A, L = triforce_refinement(EP, V, T, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = resolution, continuous = false)
            end
        end
	end
    if continuous == true
        myplot = draw_triangular_mesh(V, T, xlims = xlims, ylims = ylims, scatter1 = scatter, edges = edges, label = label, continuous = true, plotting_standard = refinement_standard)
    else
        myplot = draw_triangular_mesh(V, T, xlims = xlims, ylims = ylims, scatter1 = scatter, edges = edges, label = label, continuous = false)
    end
	return myplot
end

#=function: Delaunay_triangulation
input:
		value_dict - Dictionary that takes parameters as keys and fibre_function output for values
ouput:
		value_dict - returns the same value_dict inputted
		Triangles - function computes a Delaunay triangulation of the points given in value_dict and returns the triangles produced
=#
function Delaunay_triangulation(value_dict)
	vertices = hcat(keys(value_dict)...)
	tri = triangulate(vertices)
	triangle_iterator = each_solid_triangle(tri)
	Triangles = []
	for T in triangle_iterator
		i,j,k = triangle_vertices(T)
		i,j,k = get_point(tri, i, j, k)
		vertex_1 = [i[1], i[2]]
		vertex_2 = [j[1], j[2]]
		vertex_3 = [k[1], k[2]]
		push!(Triangles, [vertex_1,vertex_2,vertex_3])
	end
	return (value_dict, Triangles)
end

#=
function: Delaunay_visualization
input:
		EP - EnumerativeProblem
		fibre_function - function on solutions to EP
		depth - number of refinement steps, used when automatic == false
		xlims - range of "x" parameter for mesh
		ylims - range of "y" parameter for mesh 
		resolution - number of parameters to be solved for during refinement; upper limit on number of triangles to be refined = floor(resolution/3)
		total_resolution - total resolution used when automatic == true
		automatic - option to run "automatic" refinement: refinement begins with a totalResolution and runs through each level of refinement using as much resolution as totalResolution permits.
		Refinement ends when totalResolution runs out. When automatic == false, the mesh is refined n = depth times with the same upper limit on number of refinement points = resolution.
		scatter - option to plot triangle vertices
		label - option to include a legend in plot
		edges - option to plot triangle edges 
		continuous - indicate whether fibre_function is discrete or continuous
		plot_proportion - in the case that fibre_function is continuous, plot_proportion is used to calculate the threshold for both refinement and plotting. For continuous fibre functions, the
		range of fibre_function values across a triangle is calculated (range across the values of the three vertices) and this range is used to determine if the triangle will be refined. For example,
		if plot_proportion is 0.7, the 70th percentile of triangle value ranges is used as the threshold for refinement; any triangle with a value range greater than the 70th percentile will be refined. 
		When refinement is complete, any triangle with a value range less than the 70th percentile will be plotted.
output:
		myplot - A plot of the mesh following all stages of refinement
=#
function Delaunay_visualization(EP::EnumerativeProblem; xlims = [-2,2], ylims = [-2,2], fibre_function = x->HomotopyContinuation.nreal(x[1]), depth = 4, resolution = 1000, total_resolution = 8*resolution, automatic = true, scatter = false, edges = false, label = true, continuous = false, plot_proportion = 0.7)
	V, T, L = initial_triangular_mesh(EP, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = resolution)
	if automatic == false
		if continuous == true
			refinement_standard = standard_for_continuous_function_refinement(V, T, plot_proportion)
			for i in 2: depth
				V, T, A, L = triforce_refinement(EP, V, T, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = resolution, continuous = true, standard_for_refinement = refinement_standard)
				V, T = Delaunay_triangulation(V)
			end
		else
			for i in 2:depth
				V, T, A, L = triforce_refinement(EP, V, T, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = resolution, continuous = false)
				V, T = Delaunay_triangulation(V)
			end
		end
	else
		total_resolution -= L
		if continuous == true
			refinement_standard = standard_for_continuous_function_refinement(V, T, plot_proportion)
			while total_resolution > 0
				V, T, A, L = triforce_refinement(EP, V, T, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = total_resolution, continuous = true, standard_for_refinement = refinement_standard)
				V, T = Delaunay_triangulation(V)
				if A == true
					break
				end
					total_resolution -= L
			end
		else
			while total_resolution > 0
				V, T, A, L = triforce_refinement(EP, V, T, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = total_resolution, continuous = false)
				V, T = Delaunay_triangulation(V)
				if A == true
					break
				end
					total_resolution -= L
			end
		end
	end
	if continuous == true
		myplot = draw_triangular_mesh(V, T, xlims = xlims, ylims = ylims, scatter1 = scatter, label = label, edges = edges, continuous = true, plotting_standard = refinement_standard)
	else
		myplot = draw_triangular_mesh(V, T, xlims = xlims, ylims = ylims, scatter1 = scatter, label = label, edges = edges, continuous = false)
	end
	return myplot
end

#=function: collect_triangles
input:
		triangulation - triangulation data structure from DelaunayTriangulation.jl
output:
		triangles - a vector containing all triangles from triangulation
=#

function collect_triangles(triangulation)
    triangle_iterator = each_solid_triangle(triangulation)
    triangles = []
    for T in triangle_iterator
        i, j, k = triangle_vertices(T)
        i, j, k = get_point(triangulation, i, j, k)
        vertex_1 = [i[1], i[2]]
        vertex_2 = [j[1], j[2]]
        vertex_3 = [k[1], k[2]]
        push!(triangles, [vertex_1, vertex_2, vertex_3])
    end
    return triangles
end

#=function: Delaunay_triangulation_with_Ruppert_refinement
input:
		EP - EnumerativeProblem
		value_dict - Dictionary that takes parameters as keys and fibre_function output for valuesDictionary that takes parameters as keys and fibre_function output for values
		resolution - limit on number of points added by refine!
		fibre_function - function on solutions to EP
output:
		value_dict - contains all points from inputted value_dict as well as the new points added by Ruppert refinement
		triangles - contains all triangles (vector of vertices) produced by Delaunay triangulation and Ruppert refinement
=#
function Delaunay_triangulation_with_Ruppert_refinement(EP::EnumerativeProblem, value_dict::Dict; resolution = 1000, fibre_function = x->HomotopyContinuation.nreal(x[1]))
	vertices_before_refinement = keys(value_dict)
	vertices = Vector{Vector{Float64}}() #need a mutable structure of vertices to be able to use refine!
	for i in vertices_before_refinement
		push!(vertices, i)
	end
	#vertices = hcat(keys(value_dict)...)
	tri = triangulate(vertices)
	number_of_vertices = num_solid_vertices(tri)
	A = get_area(tri)
	refine!(tri, min_angle = 27, max_area = A/number_of_vertices, max_points = (number_of_vertices + resolution))
	vertex_iterator = each_solid_vertex(tri)
	vertices_added_by_refinement = []
	for V in vertex_iterator
		x = get_point(tri, V)
		this_vertex = [x[1], x[2]]
		if (this_vertex in vertices_before_refinement) == false
			push!(vertices_added_by_refinement, this_vertex)
		else
			continue
		end
	end
	println("Number of vertices added by Ruppert refinement: ", length(vertices_added_by_refinement))
	S = solve_over_params(EP, vertices_added_by_refinement, checks = [])
	for s in S
        if length(solutions(s[1]))!=degree(EP) || nsingular(s[1])!=0
        	value_dict[s[2]] = false #parameters are given false value if they produce an "error"
		else
			value_dict[s[2]]=fibre_function(s)
		end
    end
	triangles = collect_triangles(tri)
	return (value_dict, triangles)
end

#=function: Delaunay_visualization_with_Ruppert_refinement
Delaunay_visualization_with_Ruppert_refinement contains two refinement steps:
1) Refines using triforce_refinement on the basis of unresolved triangles to improve visualization of the discriminant
2) Refines using Ruppert refinement from DelaunayTriangulation.jl to improve the quality of the triangulation
input:
		EP - EnumerativeProblem
		fibre_function - function on solutions to EP
		depth - number of refinement steps, used when automatic == false
		xlims - range of "x" parameter for mesh
		ylims - range of "y" parameter for mesh 
		resolution - number of parameters to be solved for during refinement; upper limit on number of triangles to be refined = floor(resolution/3)
		total_resolution - total resolution used when automatic == true
		automatic - option to run "automatic" refinement: refinement begins with a totalResolution and runs through each level of refinement using as much resolution as totalResolution permits.
		Refinement ends when totalResolution runs out. When automatic == false, the mesh is refined n = depth times with the same upper limit on number of refinement points = resolution.
		scatter - option to plot triangle vertices
		label - option to include a legend in plot
		edges - option to plot triangle edges 
		continuous - indicate whether fibre_function is discrete or continuous
		plot_proportion - in the case that fibre_function is continuous, plot_proportion is used to calculate the threshold for both refinement and plotting. For continuous fibre functions, the
		range of fibre_function values across a triangle is calculated (range across the values of the three vertices) and this range is used to determine if the triangle will be refined. For example,
		if plot_proportion is 0.7, the 70th percentile of triangle value ranges is used as the threshold for refinement; any triangle with a value range greater than the 70th percentile will be refined. 
		When refinement is complete, any triangle with a value range less than the 70th percentile will be plotted.
output:
		myplot - A plot of the mesh following all stages of refinement
=#
function Delaunay_visualization_with_Ruppert_refinement(EP::EnumerativeProblem; xlims = [-2,2], ylims = [-2,2], fibre_function = x->HomotopyContinuation.nreal(x[1]), depth = 4, resolution = 1000, total_resolution = 8*resolution, automatic = true, scatter = false, edges = false, label = true, continuous = false, plot_proportion = 0.7)
	V, T, L = initial_triangular_mesh(EP, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = resolution)
	if automatic == false
		if continuous == true
			refinement_standard = standard_for_continuous_function_refinement(V, T, plot_proportion)
			for i in 2:depth
				V, T, A, L = triforce_refinement(EP, V, T, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = resolution, continuous = true, standard_for_refinement = refinement_standard)
				V, T = Delaunay_triangulation_with_Ruppert_refinement(EP, V, resolution = resolution, fibre_function = fibre_function)
			end
		else
			for i in 2:depth
				V, T, A, L = triforce_refinement(EP, V, T, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = resolution, continuous = false)
				V, T = Delaunay_triangulation_with_Ruppert_refinement(EP, V, resolution = resolution, fibre_function = fibre_function)
			end
		end
	else
		total_resolution -= L
		if continuous == true
			refinement_standard = standard_for_continuous_function_refinement(V, T, plot_proportion)
			while total_resolution > 0
				V, T, A, L = triforce_refinement(EP, V, T, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = total_resolution, continuous = true, standard_for_refinement = refinement_standard)
				V, T = Delaunay_triangulation_with_Ruppert_refinement(EP, V, resolution = resolution, fibre_function = fibre_function)
				if A == true
					break
				end
					total_resolution -= L
			end
		else
			while total_resolution > 0
				V, T, A, L = triforce_refinement(EP, V, T, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = total_resolution, continuous = false)
				V, T = Delaunay_triangulation_with_Ruppert_refinement(EP, V, resolution = resolution, fibre_function = fibre_function)
				if A == true
					break
				end
					total_resolution -= L
			end
		end
	end
	if continuous == true
		myplot = draw_triangular_mesh(V, T, xlims = xlims, ylims = ylims, scatter1 = scatter, label = label, edges = edges, continuous = true, plotting_standard = refinement_standard)
	else
		myplot = draw_triangular_mesh(V, T, xlims = xlims, ylims = ylims, scatter1 = scatter, label = label, edges = edges, continuous = false)
	end
	return myplot
end

#=function: standard_for_continuous_function_refinement
input:
		value_dict - Dictionary that takes parameters as keys and fibre_function output for values
		triangles - Vector containing all triangles in the triangulation
		plot_proportion - Desired proportion for determining the plotting threshold. E.g. if plot_proportion is 0.7, the 70th percentile of 
		triangle value ranges is used as the threshold for determining which triangles must be refined/plotted.
output:
		standard_for_refinement - the corresponding percentile of triangle value ranges that is used as a refinement/plotting threshold
=#
function standard_for_continuous_function_refinement(value_dict::Dict, triangles::Vector, plot_proportion = 0.7)
    triangle_value_ranges = Dict()
    for T in triangles
        min_value = min([value_dict[x] for x in T]...)
        max_value = max([value_dict[x] for x in T]...)
        triangle_range = max_value - min_value
        triangle_value_ranges[T] = triangle_range
    end
    standard_for_refinement = Statistics.quantile(values(triangle_value_ranges), plot_proportion)
    return standard_for_refinement
end