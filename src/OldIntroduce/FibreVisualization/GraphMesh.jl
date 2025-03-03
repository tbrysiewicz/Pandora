#cache's parameter->f(parameter) input/output pairs for some fibre function
#  The vector function_cache should NEVER be reordered, as Subdivision refers
#  to its elements via their indices
mutable struct GraphMesh
    function_cache :: Vector{Tuple{Vector{Float64},Float64}}
end

function input_points(GM::GraphMesh)
    [f[1] for f in GM.function_cache]
end

function output_values(GM::GraphMesh)
    [f[2] for f in GM.function_cache]
end

#TODO: Slightly misleading name since it returns the pair and not just hte value
function value_from_index(GM_index, GM::GraphMesh) #takes an index value from a subdivision polygon and returns parameter and fibre function value
    function_cache = GM.function_cache
    (parameter, parameter_fibre_function_value) = function_cache[GM_index]
    return (parameter, parameter_fibre_function_value)
end



#TODO: Find a better name. Subdivision is part of this information, but the cached
#      values make it more than just a subdivision. Maybe a ValuedSubdivision?

mutable struct Subdivision
    GM :: GraphMesh
    Polygons :: Vector{Vector{Int64}} #Indices 
end

function graph_mesh(SD::Subdivision)
    return(SD.GM)
end

function polygons(SD::Subdivision)
    return(SD.Polygons)
end