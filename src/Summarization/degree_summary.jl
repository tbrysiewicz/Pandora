function degree_summary(EP::EnumerativeProblem;kwargs...)
    summary = raw"\section{Degree Bounds, Newton Polytopes, and Support}"
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
    # Add support visualization if number of variables is 2
    if ambient_dimension(EP) == 2
        for (i, fig) in enumerate(visualize_support(EP))
            save(fig, "Newton_Polytopes$(i)")
        end
        summary *= raw"""
The Newton polytopes and support of the defining polynomials of the enumerative problem are shown in Figure~\ref{fig:NewtonPolytopes}."""

    summary *= raw"""With respect to its support, the enumerative problem is """
    if is_lacunary(EP) && is_triangular(EP)
        summary *= raw"triangular and lacunary."
    elseif is_lacunary(EP) &&  !is_triangular(EP)
        summary *= raw"lacunary but not triangular."
    elseif !is_lacunary(EP) && is_triangular(EP)
        summary *= raw"triangular but not lacunary."
    else
        summary *= raw"neither triangular nor lacunary."
    end
    summary *= raw"""

\begin{figure}[!htpb]
\centering
\includegraphics[scale=0.3]{OutputFiles/Newton_Polytopes1.png}
\includegraphics[scale=0.3]{OutputFiles/Newton_Polytopes2.png}
\caption{The Newton polytopes and support of the defining polynomials of the enumerative problem.}
\label{fig:NewtonPolytopes}
\end{figure}
"""
    end


    return(summary)
end
