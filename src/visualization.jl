export
    visualize_support,
    visualize_newton_polytopes

function vertices(P::Polyhedron)
    Vector{Vector{Int}}(1.0 * collect(oscar_vertices(PointVector, P)))
end

function dim(P::Polyhedron)
    oscar_dim(P)
end

function ambient_dimension(P::Polyhedron)
    oscar_ambient_dim(P)
end

function visualize_newton_polytopes(EP::EnumerativeProblem)
    # Compute the Newton polytopes of the enumerative problem
    NP = newton_polytopes(EP)
    # Visualize each Newton polytope
    plots = [visualize(P) for P in NP]
    return plots
end

function visualize(Ps::Vector{Polyhedron})
    # Check if all polyhedra in the vector have the same dimension
    if all(ambient_dimension(P) == ambient_dimension(Ps[1]) for P in Ps)
        #Create a plot for each polyhedron in the vector
        (Xmin,Xmax,Ymin,Ymax) = lims(Ps[1])
        for P in Ps
            (xm,xM,ym,yM) = lims(P)
            Xmin = min(Xmin, xm)
            Xmax = max(Xmax, xM)
            Ymin = min(Ymin, ym)
            Ymax = max(Ymax, yM)
        end
        plots = [visualize(P; xlims = [Xmin,Xmax], ylims = [Ymin,Ymax],aspect_ratio = :equal) for P in Ps]
        return(plots)
    else
        throw(ArgumentError("All polyhedra must have the same dimension for visualization."))
    end
end

function visualize(P::Polyhedron; kwargs...)
    #Check that dim(P) is 2 or 3, and route to the appropriate function
    #Otherwise, throw an error
    if ambient_dimension(P) == 2
        return visualize_2d_polyhedron(P; kwargs...)
    elseif ambient_dimension(P) == 3
        throw(NotImplementedError("Visualization for 3D polyhedra is not yet implemented."))
    end
    error("Visualization is only supported for 2D and 3D polyhedra.")
end

function order_vertices(edges)
    V = unique(vcat(edges...))
    ordered_vertices = edges[1]
    while length(ordered_vertices) < length(V)
        next_edge = edges[findfirst(e -> in(last(ordered_vertices), e) && !in(ordered_vertices[end-1],e),edges)]

        next_vertex = next_edge[findfirst(x->x!=last(ordered_vertices), next_edge)]
        push!(ordered_vertices, next_vertex)
    end
   return ordered_vertices
end

function lattice_plot!(xmin, xmax, ymin, ymax)
    # Initialize plot as a scatter of the relevant integer lattice, with no axes and no grid
    points = collect(Iterators.product(xmin:xmax, ymin:ymax))
    X = getindex.(points, 1)
    Y = getindex.(points, 2)    
    MyPlot = scatter!(X, Y; markersize = 3, color = :gray, legend = false, grid = false, axis = false)
    return MyPlot
end

function lims(P::Polyhedron)
    V = vertices(P)

    xmin = min(getindex.(V, 1)...)
    xmax = max(getindex.(V, 1)...)
    ymin = min(getindex.(V, 2)...)
    ymax = max(getindex.(V, 2)...)

    return(xmin, xmax, ymin, ymax)
end

function visualize_2d_polyhedron(P::Polyhedron; xlims = nothing, ylims = nothing)
    V = vertices(P)
    edges = map(vertices,faces(P,1))
    ordered_vertices = order_vertices(edges)
    
    (xmin, xmax, ymin, ymax) = lims(P)
    if xlims != nothing
        xmin, xmax = xlims
    end
    if ylims != nothing
        ymin, ymax = ylims
    end
    # Draw thick black arrows from the origin to the axes
    arrow_length = 0.2  # Adjust arrowhead size if needed
    plot([0, 0], [0, ymax+1], arrow = (:closed, 0.5), color = :black, linewidth = 3)
    plot!([0, xmax+1], [0, 0], arrow = (:closed, 0.5), color = :black, linewidth = 3)
    plot!([0, 0], [0, -1], arrow = (:closed, 0.5), color = :black, linewidth = 3)
    plot!([0, -1], [0, 0], arrow = (:closed, 0.5), color = :black, linewidth = 3)

    MyPlot = plot!(Shape(getindex.(ordered_vertices,1), getindex.(ordered_vertices,2)), 
         fillcolor = :cyan, linecolor = :blue, linewidth = 2, seriestype = :shape, aspect_ratio = :equal)


 
    lattice_plot!(-1, xmax+1, -1, ymax+1)

    return MyPlot
end


function visualize_support(EP::EnumerativeProblem)
    NP = newton_polytopes(EP)
    V = map(vertices, NP)
    NP_plots = visualize(newton_polytopes(EP))
    S = support(EP)
    S = map(s->collect(eachrow(s)), S)

    #S = [filter(s->!in(s,V[i]),S[i]) for i in eachindex(S)]
    for i in eachindex(S)
        s = S[i]
        scatter!(NP_plots[i], getindex.(s,1), getindex.(s,2); 
            markersize = 10, color = :blue)
    end
    return NP_plots
end

