import Base: getindex, iterate

using Plots
using DelaunayTriangulation: triangulate, each_solid_triangle, triangle_vertices, get_point, convert_boundary_points_to_indices
using LinearAlgebra: qr

#==============================================================================#
# EXPORTS
#==============================================================================#

export
    visualize,                   # Main user function for plotting
    visualize_function_cache,    # Plot from a function cache
    ValuedSubdivision,           # Main struct for subdivision/visualization
    refine!,                     # Refine the subdivision
    delaunay_retriangulate!,     # Re-triangulate the mesh
    is_discrete,                 # Heuristic check if function cache represents a discrete valued function
    is_complete,                 # Check if a polygon is complete
    VISUALIZATION_STRATEGIES,    # Dictionary of visualization strategies
    complete_polygons,
    incomplete_polygons,
    animate_refinement         # Animate the refinement process

#==============================================================================#
# MESH FUNCTIONS
#==============================================================================#

function initial_parameter_distribution(; kwargs...)
    xlims = get(kwargs, :xlims, [-1,1])
    ylims = get(kwargs, :ylims, [-1,1])
    resolution = get(kwargs, :resolution, 1000)
    xlength = xlims[2] - xlims[1]
    ylength = ylims[2] - ylims[1]

    x_values = range(xlims[1], xlims[2], Int(floor(sqrt((((xlength)/(xlength + ylength))*resolution)/(1 - ((xlength)/(xlength + ylength)))))))
    y_values = range(ylims[1], ylims[2], Int(floor(sqrt((((ylength)/(xlength + ylength))*resolution)/(1 - ((ylength)/(xlength + ylength)))))))

    return x_values, y_values
end

function trihexagonal_mesh(; kwargs...)
    xlims = get(kwargs, :xlims, [-1,1])
    ylims = get(kwargs, :ylims, [-1,1])
    resolution = get(kwargs, :resolution, 1000)

    x_values, y_values = initial_parameter_distribution(; xlims=xlims, ylims=ylims, resolution=resolution)

    shift_amount = (x_values[2] - x_values[1])/2 	
    parameters = [[i - isodd(findfirst(x->x==j, y_values))*shift_amount, j] for i in x_values for j in y_values]

    parameters_as_matrix = reshape(parameters, length(x_values), length(y_values))

    triangles = Vector{Vector{Int}}([])
    for i in 1:length(x_values)-1
        for j in 1:length(y_values)-1
            tl_index = findfirst(x->x == parameters_as_matrix[i,j], parameters)
            tr_index = findfirst(x->x == parameters_as_matrix[i+1,j], parameters)
            bl_index = findfirst(x->x == parameters_as_matrix[i, j+1], parameters)
            br_index = findfirst(x->x == parameters_as_matrix[i+1, j+1], parameters)

            if isodd(i)
                push!(triangles,[tl_index, tr_index, bl_index])
                push!(triangles,[bl_index, tr_index, br_index])
            else
                push!(triangles,[tl_index, tr_index, br_index])
                push!(triangles,[bl_index, tl_index, br_index])
            end
        end
    end
    return(triangles,parameters)
end

function rectangular_mesh(; kwargs...)
    xlims = get(kwargs, :xlims, [-1,1])
    ylims = get(kwargs, :ylims, [-1,1])
    resolution = get(kwargs, :resolution, 1000)
    x_values, y_values = initial_parameter_distribution(; xlims=xlims, ylims=ylims, resolution=resolution)
    parameters = [[i,j] for i in x_values for j in y_values]
    parameters_as_matrix = reshape(parameters, length(x_values), length(y_values))
    rectangles = Vector{Vector{Int}}([])
    for i in 1:length(x_values)-1
        for j in 1:length(y_values)-1
            tl_index = findfirst(x->x == parameters_as_matrix[i,j], parameters)
            tr_index = findfirst(x->x == parameters_as_matrix[i+1,j], parameters)
            bl_index = findfirst(x->x == parameters_as_matrix[i, j+1], parameters)
            br_index = findfirst(x->x == parameters_as_matrix[i+1, j+1], parameters)
            push!(rectangles, [tl_index, tr_index, br_index, bl_index])
        end
    end
    return (rectangles, parameters)
end

#==============================================================================#
# LOCAL INSERTION FUNCTIONS
#==============================================================================#

function quadtree_insertion(P::Vector{Vector{Float64}})
    center = midpoint(P)
    midpoints = [midpoint([P[i], P[mod1(i+1, 4)]]) for i in 1:4]
    rectangles = [[P[i], midpoints[i], center, midpoints[mod1(i-1, 4)]] for i in 1:4]
    return (midpoints..., center), rectangles
end

function random_point_insertion(P::Vector{Vector{Float64}})
    c = abs.(randn(Float64, length(P)))
    c ./= sum(c)
    random_point = sum(P[i] .* c[i] for i in eachindex(P))
    triangles = [[P[1], P[2], random_point],
                 [P[2], P[3], random_point],
                 [P[3], P[1], random_point]]
    return [random_point], triangles
end

function sierpinski_point_insertion(P::Vector{Vector{Float64}})
    m1 = midpoint([P[1], P[2]])
    m2 = midpoint([P[1], P[3]])
    m3 = midpoint([P[2], P[3]])
    triangles = [[P[1], m1, m2],
                 [P[2], m3, m1],
                 [P[3], m2, m3],
                 [m1, m2, m3]]
    return [m1, m2, m3], triangles
end

function barycenter_point_insertion(P::Vector{Vector{Float64}})
    barycenter = midpoint(P)
    triangles = [[P[1], P[2], barycenter],
                 [P[2], P[3], barycenter],
                 [P[3], P[1], barycenter]]
    return [barycenter], triangles
end

#==============================================================================#
# VISUALIZATION STRATEGIES
#==============================================================================#

"""
    VISUALIZATION_STRATEGIES

A dictionary of visualization strategies, each containing a mesh function and a refinement method.
"""
global VISUALIZATION_STRATEGIES = Dict{Symbol, Dict{Symbol, Any}}(
    :quadtree => Dict{Symbol, Any}(:mesh_function => rectangular_mesh, :refinement_method => quadtree_insertion),
    :barycentric => Dict{Symbol, Any}(:mesh_function => trihexagonal_mesh, :refinement_method => barycenter_point_insertion),
    :sierpinski => Dict{Symbol, Any}(:mesh_function => trihexagonal_mesh, :refinement_method => sierpinski_point_insertion),
    :random => Dict{Symbol, Any}(:mesh_function => trihexagonal_mesh, :refinement_method => random_point_insertion),
    :careful => Dict{Symbol, Any}(:mesh_function => trihexagonal_mesh, :refinement_method => barycenter_point_insertion)
)

#==============================================================================#
# FUNCTIONS NEEDED TO INITIALIZE VALUEDSUBDIVISION
#==============================================================================#

"""
    is_discrete(function_cache::Vector{Tuple{Vector{Float64},Any}})

Heuristically check if the output values in the function cache are discrete. A function is considered discrete if it has fewer than 50 unique output values.
"""
function is_discrete(function_cache::Vector{Tuple{Vector{Float64},Any}})
    # Check if the output values are discrete
    output_values = getindex.(function_cache, 2)
    unique_count = length(unique(output_values))
    return unique_count < 50
end



"""
    is_complete(p::Vector{Int}, FC::Vector{Tuple{Vector{Float64},Any}}; tol=0.0, kwargs...) 

    The default function to determine if a polygon is complete.
    A polygon is considered complete if all its vertices have been evaluated and the output values are within a specified tolerance.
"""
function is_complete(polygon::Vector{Int}, FC::Vector{Tuple{Vector{Float64},Any}}; tol = 0.0) 
    vals = [FC[v][2] for v in polygon]
    vertex_function_values = sort(filter(x->isa(x,Number),vals))
    if length(vertex_function_values) == 0
        return false # If there are no vertices, we consider it incomplete
    end
    if (vertex_function_values[end] - vertex_function_values[1]) <= tol
        return true
    else
        return false
    end
end

function default_is_complete(function_cache::Vector{Tuple{Vector{Float64},Any}})
    # Checking whether the function is continuous or discrete

    is_disc = is_discrete(function_cache)
    tol = 0.0
    if !is_disc
        values = filter(x -> isa(x, Number), getindex.(function_cache, 2))
        tol = (max(values...) - min(values...))/16
    end
    local_is_complete = (p::Vector{Int}, FC:: Vector{Tuple{Vector{Float64},Any}};kwargs...) -> is_complete(p, FC; tol=tol, kwargs...) # Default function to determine completeness of polygons
    return local_is_complete
end

#==============================================================================#
# VALUEDSUBDIVISION STRUCTURE
#==============================================================================#

"""
    ValuedSubdivision

    A mutable struct representing a subdivision of a function's domain into polygons and a caching of the function values over the vertices of 
        those polygons. 
    The fields of a ValuedSubdivision include:
        - `function_oracle`: A function that takes many parameters and returns a vector of values.
        - `function_cache`: A vector of tuples, each containing a vector of input parameters and their corresponding output value.
        - `complete_polygons`: A vector of vectors, each containing indices of the function_cache that represent complete polygons.
        - `incomplete_polygons`: A vector of vectors, each containing indices of the function_cache that represent incomplete polygons.
        - `is_complete`: A function that checks if a polygon is complete,
"""
mutable struct ValuedSubdivision
    function_oracle::Function                           #This takes MANY parameters and returns a vector of values
    function_cache::Vector{Tuple{Vector{Float64},Any}}
    complete_polygons::Vector{Vector{Int64}}            # Given as indices of function_cache
    incomplete_polygons::Vector{Vector{Int64}}
    is_complete::Function

    function ValuedSubdivision(function_oracle::Function; kwargs...)
        #Pull kwargs
        xlims = get(kwargs, :xlims, [-1,1])
        ylims = get(kwargs, :ylims, [-1,1])
        resolution = get(kwargs, :resolution, 1000)
        strategy = get(kwargs, :strategy, :sierpinski)
        mesh_function = get(kwargs, :mesh_function, VISUALIZATION_STRATEGIES[strategy][:mesh_function])

    
        VSD = new()
        VSD.function_oracle = function_oracle

        # Generate initial mesh
        (polygons, parameters) = mesh_function(; xlims=xlims, ylims=ylims, resolution=resolution, kwargs...)

        # Solve for all pts in the initial mesh
        function_oracle_values = function_oracle(parameters)
        length(function_oracle_values) != length(parameters) && error("Did not solve for each parameter")

        # Create function cache
        function_cache = Vector{Tuple{Vector{Float64},Any}}([])
        for i in eachindex(function_oracle_values)
            push!(function_cache, (parameters[i], function_oracle_values[i]))
        end
        VSD.function_cache = function_cache
        VSD.complete_polygons = Vector{Vector{Int64}}([])
        VSD.incomplete_polygons = Vector{Vector{Int64}}(polygons)

        # If is_complete was not passed, determine it. 

        if haskey(kwargs, :is_complete)
            ic = kwargs[:is_complete]
            VSD.is_complete =  (p::Vector{Int}, FC:: Vector{Tuple{Vector{Float64},Any}};kwargs...) -> ic(p, FC;kwargs...) # If a custom is_complete function is provided, use it
        else
            VSD.is_complete = default_is_complete(function_cache)
        end
        
        check_completeness!(VSD)
        return VSD
    end

    function ValuedSubdivision(EP::EnumerativeProblem; kwargs...)
        #Check enumerative problem is visualizable

        if n_parameters(EP) > 2
            near = real(get(kwargs, :near, rand(Float64,n_parameters(EP)))) # If near is not provided, we sample a random point in the parameter space
            println("EP consists of more than two parameters. Visualizing a random 2-plane in the parameter space.")
            new_EP = restrict(EP, [near, near + rand(Float64,length(near)), near + rand(Float64,length(near))])
            return(ValuedSubdivision(new_EP; kwargs...))
        end
        fibre_function = get(kwargs, :fibre_function, x->valid_real_solution_set(EP,x) ? n_real_solutions(x) : :wild)
        function_oracle = get(kwargs, :function_oracle, x->fibre_function.(EP(x)))

        ValuedSubdivision(function_oracle; kwargs...)
    end
end


function Base.show(io::IO, VSD::ValuedSubdivision)
    print(io, "ValuedSubdivision with ", length(VSD.function_cache), " function cache entries, ",
          length(VSD.complete_polygons), " complete polygons, and ",
          length(VSD.incomplete_polygons), " incomplete polygons.")
end

#==============================================================================#
# GETTERS AND SETTERS
#==============================================================================#

function_oracle(VSD::ValuedSubdivision) = VSD.function_oracle
function_cache(VSD::ValuedSubdivision) = VSD.function_cache
"""
    complete_polygons(VSD::ValuedSubdivision)

    A vector of vectors, each containing indices of the function_cache that represent complete polygons.
"""
complete_polygons(VSD::ValuedSubdivision) = VSD.complete_polygons

"""
    incomplete_polygons(VSD::ValuedSubdivision)

    A vector of vectors, each containing indices of the function_cache that represent incomplete polygons.
"""
incomplete_polygons(VSD::ValuedSubdivision) = VSD.incomplete_polygons
input_points(FC::Vector{Tuple{Vector{Float64}, Any}}) = getindex.(FC, 1)
output_values(FC::Vector{Tuple{Vector{Float64}, Any}}) = getindex.(FC, 2)
input_points(VSD::ValuedSubdivision) = input_points(function_cache(VSD))
output_values(VSD::ValuedSubdivision) = output_values(function_cache(VSD))
is_discrete(VSD::ValuedSubdivision) = is_discrete(function_cache(VSD))

function is_complete(polygon::Vector{Int}, VSD::ValuedSubdivision; kwargs...) 
    VSD.is_complete(polygon, function_cache(VSD); kwargs...)
end

#==============================================================================#
# POLYGON MANAGEMENT
#==============================================================================#

function delete_from_incomplete_polygons!(VSD::ValuedSubdivision, P::Vector{Int64}; kwargs...)
    polygon_index = findfirst(x -> x == P, incomplete_polygons(VSD))
    if polygon_index !== nothing
        deleteat!(VSD.incomplete_polygons, polygon_index)
    else
        error("Polygon not found in IncompletePolygons!")
    end
end

function push_to_complete_polygons!(VSD::ValuedSubdivision, P::Vector{Int64}; kwargs...)
    push!(VSD.complete_polygons, P)
end

function check_completeness!(VSD::ValuedSubdivision; kwargs...)
    # Check if the current polygons are complete
    complete_bucket = []
    for p in incomplete_polygons(VSD)
        if is_complete(p, VSD; kwargs...)
            push!(complete_bucket, p)
        end
    end
    for p in complete_bucket
        delete_from_incomplete_polygons!(VSD, p)
        push_to_complete_polygons!(VSD, p)
    end
end

function set_complete_polygons!(VSD::ValuedSubdivision, P::Vector{Vector{Int64}})
    VSD.complete_polygons = P
end

function set_incomplete_polygons!(VSD::ValuedSubdivision, P::Vector{Vector{Int64}})
    VSD.incomplete_polygons = P
end

function push_to_incomplete_polygons!(VSD::ValuedSubdivision, P::Vector{Int64}; kwargs...)
    push!(VSD.incomplete_polygons, P)
end

function push_to_function_cache!(VSD::ValuedSubdivision, V::Tuple{Vector{Float64},Any}; kwargs...)
    push!(VSD.function_cache, V)
end

#==============================================================================#
# REFINEMENT FUNCTION
#==============================================================================#

function area_of_polygon(P::Vector{Vector{Float64}})::Float64
    # Computes the area of a simple polygon using the shoelace formula
    n = length(P)
    area = 0.0
    for i in 1:n
        x1, y1 = P[i]
        x2, y2 = P[mod1(i+1, n)]
        area += x1 * y2 - x2 * y1
    end
    return 0.5 * abs(area)
end

function refine!(VSD::ValuedSubdivision, resolution::Int64;	strategy = nothing, kwargs...)

    if strategy === nothing
        n = length(complete_polygons(VSD)[1])
        strategy = n == 4 ? :quadtree : n == 3 ? :sierpinski : strategy
    end

    if in(strategy, collect(keys(VISUALIZATION_STRATEGIES))) == false
        error("Invalid strategy inputted. Valid strategies include:",keys(VISUALIZATION_STRATEGIES),".")
    end

    local_refinement_method = get(kwargs, :refinement_method, VISUALIZATION_STRATEGIES[strategy][:refinement_method])

    FO = function_oracle(VSD)
    refined_polygons::Vector{Vector{Int64}} = []
    polygons_to_solve_and_sort::Vector{Vector{Vector{Float64}}} = []
    new_parameters_to_solve::Vector{Vector{Float64}} = []
    resolution_used = 0


    IP = sort(incomplete_polygons(VSD), by = x -> area_of_polygon([VSD.function_cache[v][1] for v in x]), rev = true)

    for P in IP
        P_parameters = [function_cache(VSD)[v][1] for v in P]
        new_params, new_polygons = local_refinement_method(P_parameters)
        push!(refined_polygons, P)
        push!(polygons_to_solve_and_sort, new_polygons...)
        for p in new_params
            if !(p in new_parameters_to_solve) && !(p in input_points(VSD))
                push!(new_parameters_to_solve, p)
                resolution_used += 1
            end
        end
        resolution_used >= resolution && break
    end

    for P in refined_polygons
        delete_from_incomplete_polygons!(VSD, P)
    end

    length(new_parameters_to_solve) == 0 && return resolution_used
    function_oracle_values = FO(new_parameters_to_solve)
    length(function_oracle_values) != length(new_parameters_to_solve) && error("Did not solve for each parameter")
    for i in eachindex(function_oracle_values)
         push_to_function_cache!(VSD, (new_parameters_to_solve[i], (function_oracle_values[i])))
    end

    for P in polygons_to_solve_and_sort
        polygon = [findfirst(x->x[1]==y, function_cache(VSD)) for y in P]
        push_to_incomplete_polygons!(VSD, polygon) 
    end
    
    check_completeness!(VSD)
    println("Resolution used:", resolution_used)
    return(VSD)
end

function refine!(VSD::ValuedSubdivision; kwargs...)
    # Default refinement strategy is quadtree
    return refine!(VSD, 1000000; kwargs...)
end

#==============================================================================#
# DELAUNAY RETRIANGULATION
#==============================================================================#

"""
    delaunay_retriangulate!(VSD::ValuedSubdivision)

    Re-triangulates the mesh of a ValuedSubdivision using Delaunay triangulation.
    This function ensures that the triangulation does not contain duplicate points and updates the polygons accordingly.
"""
function delaunay_retriangulate!(VSD::ValuedSubdivision)
    @vprintln("Delaunay re-triangulating the mesh...")
    vertices = []
    for v in function_cache(VSD)
        v[1] in vertices || push!(vertices, v[1]) #ensuring that the triangulation will not contain duplicate points and the package won't give that annoying warning
    end
    vertices = hcat(vertices...)
    tri = triangulate(vertices) #This comes from DelaunayTriangulation.jl
    triangle_iterator = each_solid_triangle(tri) #Also from DelaunayTriangulation.jl
    triangles::Vector{Vector{Int64}} = []
    for T in triangle_iterator
        i,j,k = triangle_vertices(T)
        i,j,k = get_point(tri, i, j, k)
        vertex_1 = [i[1], i[2]]
        vertex_2 = [j[1], j[2]]
        vertex_3 = [k[1], k[2]]
        index_1 = findfirst(x->x[1] == vertex_1, function_cache(VSD))
        index_2 = findfirst(x->x[1] == vertex_2, function_cache(VSD))
        index_3 = findfirst(x->x[1] == vertex_3, function_cache(VSD))
        triangle = [index_1,index_2,index_3]
        push!(triangles, triangle)
    end
    set_complete_polygons!(VSD, Vector{Vector{Int64}}([])) 
    set_incomplete_polygons!(VSD, triangles)
    check_completeness!(VSD)
    return VSD
end

#==============================================================================#
# UTILITY FUNCTIONS
#==============================================================================#

function midpoint(P::AbstractVector) 
    if length(P) == 0
        return nothing
    end
    return((sum(P))./length(P))
end

mean(P::Vector{T}) where T <: Any = midpoint(P)

function median(V::AbstractVector)
    Vsorted = sort(V)
    n = length(Vsorted)
    if n == 0
        error("Cannot compute median of empty vector")
    elseif isodd(n)
        return Vsorted[(n+1) ÷ 2]
    else
        return (Vsorted[n ÷ 2] + Vsorted[n ÷ 2 + 1]) / 2
    end
end

#==============================================================================#
# VISUALIZATION FUNCTIONS
#==============================================================================#

"""
    visualize_function_cache(VSD::ValuedSubdivision)

    Visualizes the function cache of a ValuedSubdivision by scatter plotting the input points and their corresponding output values.
    This function is useful for quickly visualizing the cached values without needing to create a full subdivision visualization.
"""
function visualize_function_cache(VSD::ValuedSubdivision)
    visualize(function_cache(VSD))
end
function visualize(FC::Vector{Tuple{Vector{Float64},Any}})
    scatter(first.(input_points(FC)), last.(input_points(FC)), 
    zcolor = map(x->isa(x,Number) ? x : -10.0, output_values(FC)), legend = false, colorbar = true)
end

"""
    visualize(VSD::ValuedSubdivision; kwargs...)
    Visualizes the ValuedSubdivision by plotting the polygons and their associated values.
    The function accepts keyword arguments for customization, such as 
        -`xlims`, 
        -`ylims`,
        -`plot_log_transform`,
        -`plot_all_polygons`
"""
function visualize(VSD::ValuedSubdivision; kwargs...)
    xl = get(kwargs, :xlims, [min(map(x->x[1][1],Pandora.function_cache(VSD))...),max(map(x->x[1][1],Pandora.function_cache(VSD))...)])
    yl = get(kwargs, :ylims, [min(map(x->x[1][2],Pandora.function_cache(VSD))...),max(map(x->x[1][2],Pandora.function_cache(VSD))...)])
    plot_log_transform = get(kwargs, :plot_log_transform, false)
    plot_all_polygons = get(kwargs, :plot_all_polygons, is_discrete(VSD) == false)

    MyPlot = plot(xlims = xl, ylims = yl, aspect_ratio = :equal, background_color_inside=:black; kwargs...)
    
    polygon_list = plot_all_polygons == true ? vcat(complete_polygons(VSD), incomplete_polygons(VSD)) : complete_polygons(VSD)
    polygons_to_draw = map(P->map(first,function_cache(VSD)[P]), polygon_list)
    # we remove non-numbers when computing means so that non-numbers function as wild cards essentially 
    polygon_values = map(P-> mean(filter(x->isa(x,Number),map(last,function_cache(VSD)[P]))), polygon_list)
    if plot_log_transform
        # If we are plotting log transformed values, we need to transform the polygon values
        polygon_values = map(x->x==nothing ? nothing : log(x+1), polygon_values)
    end
    # It is still possible that the mean returns 'nothing' if all values are non-numbers
    real_values = filter(x->x!=nothing, unique(polygon_values))
    max_value = max(real_values...)
    shapes = Shape[]
    colors = []

    for (poly_pts, value) in zip(polygons_to_draw, polygon_values)
        push!(shapes, Shape([(t[1], t[2]) for t in poly_pts]))
        push!(colors, value==nothing ? :white : cgrad(:thermal, rev=false)[value/max_value])
    end
    MyPlot = plot!(shapes; fillcolor=permutedims(colors), linecolor=permutedims(colors), linewidth=0, label = false)

    
    legend_values = sort(real_values, rev=true)
    if in(nothing,polygon_values)
        pushfirst!(legend_values, nothing)
    end
    if length(legend_values) <= 10
        # Add legend entries for unique values
        for val in legend_values
            color = val==nothing ? :white : cgrad(:thermal, rev=false)[val/max_value]
            # Plot an invisible shape with the correct color and label for the legend
            plot!(Shape([(NaN, NaN), (NaN, NaN), (NaN, NaN)]); fillcolor=color, linecolor=:transparent, linewidth=0, label="$val", legend = :outerright)
        end
    else
        # Add a colorbar for the continuous values
        scatter!([NaN], [NaN]; zcolor = [min(real_values...), max(real_values...)], color = cgrad(:thermal, rev = false), markersize = 0, colorbar = true, label = false)
    end
    return(MyPlot)
end

"""
    visualize(EP::EnumerativeProblem; kwargs...)
    visualize(EP::EnumerativeProblem; fibre_function = ..., kwargs...)

    Visualizes an EnumerativeProblem `EP` by plotting the values of `fibre_function`, applied to fibres of `EP`, over the parameter space.
    If the parameter space has more than two parameters, it restricts the visualization to a random 2-plane in the parameter space.
    
    The default fibre function is `n_real_solutions`, which counts the number of real solutions in the fibre.
    Other interesting fibre functions include 
    - `dietmaier`

    Useful keyword arguments include:
    - `xlims`: Limits for the x-axis.
    - `ylims`: Limits for the y-axis.
    - `resolution`: Number of points to sample in the parameter space.
    - `strategy`: Strategy for mesh generation, e.g., `:quadtree`, `:barycentric`, `:sierpinski`, or `:random`.
    - `plot_log_transform`: Whether to apply a logarithmic transformation to the output values for better visualization.
    - `plot_all_polygons`: Whether to plot all polygons or only the complete ones. 
    - `near`: A vector of parameters to visualize around, useful for zooming in on a specific region of the parameter space.
"""
function visualize(EP::EnumerativeProblem; kwargs...)

    #Check enumerative problem is visualizable
    if n_parameters(EP) > 2
        near = real(get(kwargs, :near, rand(Float64,n_parameters(EP)))) # If near is not provided, we sample a random point in the parameter space
        println("EP consists of more than two parameters. Visualizing a random 2-plane in the parameter space.")
        new_EP = restrict(EP, [near, near + rand(Float64,length(near)), near + rand(Float64,length(near))])
    else
        new_EP = EP
    end

    strat = get(kwargs, :strategy, :careful)

    VSD = ValuedSubdivision(new_EP; kwargs...)

    repeats = strat == :careful ? 2 : 2
    for i in 1:repeats
        refine!(VSD;kwargs...)
        if strat == :careful
            delaunay_retriangulate!(VSD)
        end
    end

    my_plot = visualize(VSD; kwargs...)
    display(my_plot)
    return VSD
end


"""
    animate_refinement(VSD::ValuedSubdivision; steps::Int=100, resolution::Int=100, kwargs...)

    Creates an animation of the refinement process of a ValuedSubdivision `VSD`.
    The animation shows the subdivision being refined step by step, with each step visualized using the `visualize` function.
    The function accepts the following keyword arguments:
    - `steps`: The number of refinement steps to animate (default is 100).
    - `resolution`: The resolution for each refinement step (default is 100).
    - `kwargs`: Additional keyword arguments to pass to the `visualize` function for customization
"""
function animate_refinement(VSD::ValuedSubdivision; steps::Int=100, resolution::Int=100, kwargs...)

    anim = @animate for i in 1:steps
        refine!(VSD, resolution)
        visualize(VSD; kwargs...)
    end

    gif(anim, "refinement_animation.gif", fps=10)
end