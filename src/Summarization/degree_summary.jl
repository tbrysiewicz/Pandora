function degree_summary(EP::EnumerativeProblem;kwargs...)
    summary = raw"\section{Degree Bounds}"
    if knows(EP,DEGREE)
        summary *= raw"The degree of the enumerative problem is $"*string(degree(EP))*raw"$. "
    end

    if knows(EP,BEZOUT_BOUND)
        summary *= raw"The B\'ezout bound of the enumerative problem is $"*string(bezout_bound(EP))*raw"$. "
    end
    
    if knows(EP,BKK_BOUND)
        summary *= raw"The BKK bound of the enumerative problem is $"*string(bkk_bound(EP))*raw"$. "
    end

    if knows(EP,AFFINE_BKK_BOUND)
        summary *= raw"The affine BKK bound of the enumerative problem is $"*string(affine_bkk_bound(EP))*raw"$. "
    end

    return(summary)
end
