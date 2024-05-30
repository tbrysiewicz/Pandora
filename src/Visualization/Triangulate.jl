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
            tl = [xrange[i],yrange[j]]
            tr = [xrange[i+1],yrange[j]]
            bl = [xrange[i],yrange[j+1]]
            br = [xrange[i+1],yrange[j+1]]
            push!(Triangles,[tl,tr,bl])
            push!(Triangles,[bl,tr,br])
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