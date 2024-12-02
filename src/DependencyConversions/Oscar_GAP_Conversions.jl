function subgroup(perms::Vector{PermGroupElem})
    H,_  = sub(symmetric_group(degree(perms[1])),perms)
    return(H)
end