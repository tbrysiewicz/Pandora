mutable struct GraphMesh
    function_cache :: Vector{Tuple{Vector{Float64},Float64}}
end

mutable struct Subdivision
    GM :: GraphMesh
    Polygons :: Vector{Vector{Int64}} #Indices 
end
