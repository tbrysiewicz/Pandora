

function monodromy_coordinates(EP::EnumerativeProblem)

end

#Should this be naive transitivity coordinates or should we implement
#  a binary tree structure to appropriately compress the information
#  (i.e. founders, cycle reps, and cycle lengths)
function transitivity_coordinates(G::Vector{PermGroupElem})
    base_point = 1

end