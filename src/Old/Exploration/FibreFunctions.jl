export
    real_points_in_fibre,
    positive_points_in_fibre




@doc raw"""
    real_points_in_fibre(E::EnumerativeProblem,PathRes::Result,Param::Vector{Vector{ComplexF64}})

 Counts the number of paths in `PathRes` which result in real solutions (i.e. every coordinate is (approximately) real). 
 """
function real_points_in_fibre(E::EnumerativeProblem,PathRes::Result,Param::Union{Vector{ComplexF64},Vector{Float64}})
    n_real(solutions(PathRes))
end

@doc raw"""
    positive_points_in_fibre(E::EnumerativeProblem,PathRes::Result,Param::Vector{Vector{ComplexF64}})

 Counts the number of paths in `PathRes` which result in positive solutions (i.e. every coordinate is (approximately) positive). 
 """
function positive_points_in_fibre(E::EnumerativeProblem,PathRes::Result,Param::Union{Vector{ComplexF64},Vector{Float64}})
    rs = HomotopyContinuation.real_solutions(PathRes)
    ps = []
    for r in rs
        positive=true
        for a in r
            if a<0
                positive=false
            end
        end
        if positive==true
            push!(ps,r)
        end
    end
    return(length(ps))
end


