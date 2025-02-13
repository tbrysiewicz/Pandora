
function is_complete(p::Vector{Int}, GM::GraphMesh) #checks if all polygon vertices share the same fibre function value. 
                                       #If polygon is complete, also returns new parameters to insert.
    v1_value = value_from_index(p[1],GM)[2]
    for i in p[2:end]
        vval = value_from_index(i,GM)[2]
        if vval!=v1_value
            return(false)
        end
    end
    return(true)
end

function midpoint(p,q)
    return((p+q)./2)
end

#TODO: how should this be done for non-triangles? Should it be done at all for non-triangles?
#TODO: this will produce duplicate midpoints on edges which are shared by multiple incomplete triangles!
function subdivide(p::Vector{Int},GM::GraphMesh)
     v = [value_from_index(p[i],GM)[1] for i in 1:3]
    return([sum(v)./3])
   # new_parameters = [midpoint(vertices[i],vertices[j]) for (i,j) in [(1,2),(1,3),(2,3)]]
   # return(new_parameters)
end

function delaunay_triangulation!(SD::Subdivision)
    vertices = []
    for v in SD.GM.function_cache
        push!(vertices, v[1])
    end
	vertices = hcat(vertices...)
	tri = triangulate(vertices)
	triangle_iterator = each_solid_triangle(tri)
	triangles = []
	for T in triangle_iterator
		i,j,k = triangle_vertices(T)
		i,j,k = get_point(tri, i, j, k)
		vertex_1 = [i[1], i[2]]
		vertex_2 = [j[1], j[2]]
		vertex_3 = [k[1], k[2]]
        index_1 = findfirst(x->x[1] == vertex_1, SD.GM.function_cache)
        index_2 = findfirst(x->x[1] == vertex_2, SD.GM.function_cache)
        index_3 = findfirst(x->x[1] == vertex_3, SD.GM.function_cache)
		push!(triangles, [index_1, index_2, index_3])
	end
    SD.Polygons = triangles
	return nothing
end

function refine!(S::Subdivision, EP::EnumerativeProblem, resolution = 1000; fibre_function = x->n_real_solutions(x))
    new_parameters = Vector{Vector{Float64}}([])
    resolution_used = 0
    for p in S.Polygons
        if resolution_used > (resolution - 3)
            break
        else
            check_polygon = is_complete(p, S.GM)
            if check_polygon == false
                some_new_parameters = subdivide(p,S.GM)
                for snp in some_new_parameters
                    push!(new_parameters,snp)
                end
                resolution_used += 3
            end
        end
    end
    new_solutions = solve(EP, new_parameters)
    for s in new_solutions
        push!(S.GM.function_cache, (new_parameters[findfirst(x->x==s, new_solutions)], fibre_function(s)))
    end

    delaunay_triangulation!(S)

    return resolution_used
end