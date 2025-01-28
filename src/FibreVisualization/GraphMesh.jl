#cache's parameter->f(parameter) input/output pairs for some fibre function
#  The vector function_cache should NEVER be reordered, as Subdivision refers
#  to its elements via their indices
mutable struct GraphMesh
    function_cache :: Vector{Tuple{Vector{Float64},Float64}}
end

mutable struct Subdivision
    GM :: GraphMesh
    Polygons :: Vector{Vector{Int64}} #Indices 
end
