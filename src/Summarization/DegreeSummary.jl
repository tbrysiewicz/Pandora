function degree_summary(EP::EnumerativeProblem;kwargs...)
    summary = raw"\section{Degree Bounds}"
    if haskey(data(EP),:degree)
        summary *= raw"The degree of the enumerative problem is $"*string(degree(EP))*raw"$. "
    end

    if haskey(data(EP),:bezout)
        summary *= raw"The B\'ezout bound of the enumerative problem is $"*string(bezout(EP))*raw"$. "
    end
    
    if haskey(data(EP),:bkk)
        summary *= raw"The BKK bound of the enumerative problem is $"*string(bkk(EP))*raw"$. "
    end

    if haskey(data(EP),:affine_bkk)
        summary *= raw"The affine BKK bound of the enumerative problem is $"*string(affine_bkk(EP))*raw"$. "
    end

end
