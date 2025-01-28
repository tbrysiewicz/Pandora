
#write function that finds function value of an index associated with a parameter from the subdivision (it needs to search through the graphmesh)

function fibre_function_value(GM_index, GM::GraphMesh) #takes an index value from a subdivision polygon and returns parameter and fibre function value
    function_cache = GM.function_cache
    (parameter, parameter_fibre_function_value) = function_cache[GM_index]
    return (parameter, parameter_fibre_function_value)
end

function is_complete(p, GM::GraphMesh) #checks if all polygon vertices share the same fibre function value. If polygon is complete, also returns new parameters to insert.
    vertex_1_fibre_function_value = fibre_function_value(p[1], GM)
    vertex_2_fibre_function_value = fibre_function_value(p[2], GM)
    vertex_3_fibre_function_value = fibre_function_value(p[3], GM)
    new_parameters = []
    if vertex_1_fibre_function_value[2] != vertex_2_fibre_function_value[2] || vertex_1_fibre_function_value[2] != vertex_3_fibre_function_value[2]
        midpoint1 = 0.5*(vertex_2_fibre_function_value[1]-vertex_1_fibre_function_value[1]) + vertex_1_fibre_function_value[1]
        midpoint2 = 0.5*(vertex_3_fibre_function_value[1]-vertex_2_fibre_function_value[1]) + vertex_2_fibre_function_value[1]
        midpoint3 = 0.5*(vertex_3_fibre_function_value[1]-vertex_1_fibre_function_value[1]) + vertex_1_fibre_function_value[1]
        push!(new_parameters, midpoint1)
        push!(new_parameters, midpoint2)
        push!(new_parameters, midpoint3)
        return (false, new_parameters)
    else
        return (true, new_parameters)
    end
end

function Delaunay_triangulation!(SD::Subdivision)
    vertices = []
    for v in SD.GM.function_cache
        push!(vertices, v[1])
    end
	vertices = hcat(vertices...)
	tri = DelaunayTriangulation.triangulate(vertices)
	triangle_iterator = each_solid_triangle(tri)
	triangles = []
	for T in triangle_iterator
		i,j,k = triangle_vertices(T)
		i,j,k = get_point(tri, i, j, k)
		vertex_1 = [i[1], i[2]]
		vertex_2 = [j[1], j[2]]
		vertex_3 = [k[1], k[2]]
        index_1 = findfirst(x->x[1] == vertex_1, SM.GM.function_cache)
        index_2 = findfirst(x->x[1] == vertex_2, SM.GM.function_cache)
        index_3 = findfirst(x->x[1] == vertex_3, SM.GM.function_cache)
		push!(triangles, [index_1, index_2, index_3])
	end
    SD.Polygons = triangles
	return nothing
end

function refine!(S::Subdivision, EP::EnumerativeProblem, resolution = 1000; fibre_function = x->n_real_solutions(x)) :: Tuple{GraphMesh,Subdivision}
    new_parameters = []
    resolution_used = 0
    for p in S.Polygons
        if resolution_used > (resolution - 3)
            break
        else
            check_polygon = is_complete(p, S.GM)
            if check_polygon[1] == false
                for i in check_polygon[2]
                    push!(new_parameters, i)
                end
                resolution_used += 3
            end
        end
    end
    new_solutions = solve(EP, new_parameters)
    for s in new_solutions
        push!(S.GM.function_cache, (new_parameters[findfirst(x->x==s, new_solutions)], fibre_function(s)))
    end

    Delaunay_triangulation!(S)

    return resolution_used
end