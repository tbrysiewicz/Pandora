

function refine!(GM::GraphMesh, S::Subdivision, EP::EnumerativeProblem, f) :: Tuple{GraphMesh,Subdivision}
    new_inputs = Vector{Vector{Float64}}()
    #Find the new parameters you need
    for polygon in polygons(S)
        if is_complete(GM,polygon)==false
            push!(new_inputs,include_more_points(GM,polygon))
        end
    end

    #Now solve for the new parameters
    
    #Now create brand new subdivision via Delaunay

    #And optionally reduce it in some clever way. 

end