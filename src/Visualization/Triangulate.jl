#using Random

using DelaunayTriangulation

export
	visualize_with_triangles,
	triforce_visualization

function visualize_with_triangles(EP;depth = 4, resolution=1000,xlims=[-2,2],ylims=[-2,2],fibre_function = x->HomotopyContinuation.nreal(x[1]), scatter = true)
    (V,T) = initial_triangular_mesh(EP;
                                    fibre_function = fibre_function, 
                                    xlims = xlims,
                                    ylims=ylims,
                                    resolution=resolution)
    for i in 2:depth
    	println("Refinement step: ", (i-1))
        (V,T) = refine_triangular_mesh(EP,V,T;
                                    fibre_function = fibre_function, 
                                    xlims = xlims,
                                    ylims=ylims,
                                    resolution=resolution)
    end

    myplot = draw_triangular_mesh(V,T, scatter1 = scatter)
end


function refine_triangular_mesh(EP,value_dict,Triangles; fibre_function = x->HomotopyContinuation.nreal(x[1]), xlims = [-2,2],ylims=[-2,2],resolution=1000)
    non_complete = filter(T->value_dict[T[1]]!=value_dict[T[2]] || value_dict[T[1]] != value_dict[T[3]],Triangles)

    #=
    if length(non_complete)>resolution
    	non_complete = unique(rand(non_complete, Int(resolution))
    	end
    		 
  
    end
	=#
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
        if length(solutions(s[1]))!=degree(EP) || nsingular(s[1])!=0
        	value_dict[s[2]] = -2
        end
    end
    return((value_dict,new_triangles))
end

function initial_triangular_mesh(EP;fibre_function = x->HomotopyContinuation.nreal(x[1]), xlims = [-2,2],ylims=[-2,2],resolution=1000)
    nxy = Int(floor(sqrt(resolution+100)))
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
        	value_dict[s[2]] = -2
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
    return((value_dict,Triangles))
end


function draw_triangle!(T,v)
    Triangle = Shape([(t[1],t[2]) for t in T])
    c =cgrad(:thermal, rev = false)[v]
    plot!(Triangle, fillcolor = c, linecolor=c, linewidth=0, labels = false)
end

function draw_triangular_mesh(value_dict,Triangles;xlims = [-2,2],ylims=[-2,2], scatter1 = true)
    M = max(values(value_dict)...)
    m = min(values(value_dict)...)
    
    myplot = plot(xlims=xlims,ylims=ylims,legend=true)
    for T in Triangles
        if value_dict[T[1]]==value_dict[T[2]] && value_dict[T[1]]==value_dict[T[3]]
            draw_triangle!(T,value_dict[T[1]]/M)
        end
    end
    
    if scatter1 == true
    	for i in unique(values(value_dict))
    		temp_parameters = filter(x->value_dict[x]==i, keys(value_dict))
    		if i == -2
    			scatter!([A[1] for A in temp_parameters], [A[2] for A in temp_parameters], markercolor =:red, markersize = 0.8, markershape =:rect, markerstrokewidth = 0.2, labels = "Error")
    		else
    		c =cgrad(:thermal, rev = false)[i/M]
    		scatter!([A[1] for A in temp_parameters], [A[2] for A in temp_parameters], markercolor = c, markersize = 0.8, markershape =:rect, markerstrokewidth = 0.2, labels = "$(i) real solutions")
    		end
    	end
    end
    		
    
    return(myplot)
end



#Below is code for "triforce" triangulation visualization strategy


function refine_triangular_mesh2(EP::EnumerativeProblem, value_dictionary, triangles; fibre_function = x->HomotopyContinuation.nreal(x[1]), xlims = [-2,2],ylims=[-2,2],resolution=1000)
	whiteSpace = filter(x->value_dictionary[x[1]]!=value_dictionary[x[2]] || value_dictionary[x[1]]!=value_dictionary[x[3]], triangles)

	newParameters = []
	
	TrianglesToRefine=[]
	if (length(whiteSpace)*3)>resolution
		while length(TrianglesToRefine)*3<resolution
			v = rand(whiteSpace,Int(floor(resolution/10)))
			TrianglesToRefine = unique(vcat(TrianglesToRefine,v))
		end
		TrianglesToRefine = TrianglesToRefine[1:Int(floor(resolution/3))]
#		whiteSpace = unique(rand(whiteSpace,Int(floor(resolution/3))))
	else
		TrianglesToRefine = whiteSpace
	end
	
	new_triangles = setdiff(triangles, TrianglesToRefine) #Triangles not to refine. 

	
		
	for i in TrianglesToRefine
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
		if length(solutions(i[1]))!=degree(EP) || nsingular(i[1])!=0
        	value_dictionary[i[2]] = -2
        	end
	end
	
	return((value_dictionary, new_triangles))
end

function refine_mesh_for_automatic(EP::EnumerativeProblem, value_dict, triangles, total_resolution; xlims = [-2,2], ylims = [-2,2], fibre_function = x->HomotopyContinuation.nreal(x[1]))
	whiteSpace = filter(x->value_dict[x[1]]!=value_dict[x[2]] || value_dict[x[1]]!=value_dict[x[3]], triangles)
	TrianglesToRefine = []
	Termination = false
	if (length(whiteSpace)*3)>total_resolution
		while length(TrianglesToRefine)*3<total_resolution
			v = rand(whiteSpace,Int(floor(total_resolution/10)))
			TrianglesToRefine = unique(vcat(TrianglesToRefine,v))
		end
		TrianglesToRefine = TrianglesToRefine[1:Int(floor(total_resolution/3))]
		Termination = true
	else
		TrianglesToRefine = whiteSpace
	end
	
	newParameters = []
	new_triangles = setdiff(triangles, TrianglesToRefine) #Triangles not to refine. 
	
	for i in TrianglesToRefine
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
		value_dict[i[2]] = fibre_function(i)
		if length(solutions(i[1]))!=degree(EP) || nsingular(i[1])!=0
        	value_dict[i[2]] = -2
        	end
	end
	
	total_resolution = total_resolution - length(newParameters)
	
	return((value_dict, new_triangles, total_resolution, Termination))
end
		

function triforce_visualization(EP;depth = 4, initial_resolution=1000, total_resolution = 10*initial_resolution, xlims=[-2,2],ylims=[-2,2],fibre_function = x->HomotopyContinuation.nreal(x[1]), scatter = true, automatic = true)
  
    if automatic == true
    	
    	(V,T) = initial_triangular_mesh(EP;
                                    fibre_function = fibre_function, 
                                    xlims = xlims,
                                    ylims=ylims,
                                    resolution=initial_resolution)
         total_resolution = total_resolution - initial_resolution
         
         while total_resolution > 0
         	(V,T,R,D) = refine_mesh_for_automatic(EP,V,T, total_resolution;
                                    fibre_function = fibre_function, 
                                    xlims = xlims,
                                    ylims=ylims)
                     
                total_resolution = R
                
		if D == true
			break
		end
	end
	
	else
   	 (V,T) = initial_triangular_mesh(EP;
                                    fibre_function = fibre_function, 
                                    xlims = xlims,
                                    ylims=ylims,
                                    resolution=initial_resolution)
   	 for i in 2:depth
    		println("Refinement step: ", (i-1))
       		 (V,T) = refine_triangular_mesh2(EP,V,T;
                                    fibre_function = fibre_function, 
                                    xlims = xlims,
                                    ylims=ylims,
                                    resolution=initial_resolution)
   	 end
   	 end

    myplot = draw_triangular_mesh(V,T, scatter1 = scatter)
end
