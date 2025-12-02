export
    reality_summary

#TODO: summarys should not run computations, but rather, report cached info
#      so explore_reality should produce some sort of knowledge structure like
#      data_sample
function reality_summary(EP::EnumerativeProblem;kwargs...)
    n_samples = get(kwargs, :n_samples, 1000)
    summary = raw"""\section{Real Solutions}"""
    sampler = UniformSampler{Float64}(n_parameters(EP))
    EX = explore_reality(EP; sampler = sampler, n_samples = n_samples, kwargs...)
    SD = SampleDatum(n_real_solutions, sampler, EX)

    possible_reals = (degree(EP) % 2):2:degree(EP)

    H = histogram(SD; bins = 0:1:degree(EP)+2, kwargs...)
    save(H, "reality_histogram.png")

    summary *= raw"""

\begin{figure}[!htpb]
\centering
\includegraphics[scale=0.3]{OutputFiles/reality_histogram.png}
\caption{Histogram of the number of real solutions found in the fibres of enumerative problem.}
\label{fig:reality_histogram}
\end{figure}
"""

    summary *= "\n\n"
    summary *= "We sampled "*string(n_samples)*" fibres, and the number of real solutions found is displayed in the histogram"*raw"""\ref{fig:reality_histogram}"""
    summary *= "\n\n"
    summary *= raw"""The tally of numbers of real solutions is """
    summary *= raw"""\begin{align*}"""*"\n"
    tally_pairs = []
    for (k,v) in tally(SD.function_values)
        if k in possible_reals
            push!(tally_pairs, (k,v))
        end
    end
    tally_pairs = sort(tally_pairs, by = x -> x[1])  # Sort by the number of real solutions
    for (k,v) in tally_pairs
        summary *= string(k)*raw"""&="""*string(v)*raw""" \\ """
    end
    summary *= raw"""\end{align*}"""*"\n"
    #take set difference with the values in EX

    O = maximize_n_real_solutions(EP;  n_samples = 20, max_steps = 10, kwargs...)
    maxO = n_real_solutions(record_fibre(O)[1])
    maxEX = max(EX...)
    max_real = max(maxO, maxEX)
    if maxO==max_real
        summary *= raw"""The maximum number of real solutions found in the fibres is $"""*string(max_real)*raw"""$. """
        summary *= raw"""This number was realized by running Pandora's optimizer on the so-called Dietmaier scheme"""
        summary *= raw""", which was run on 10 steps of 20 samples. It is realized by the fibre"""
        summary *= raw""" \begin{align*}"""*"\n"
        for i in 1:n_parameters(EP)
            p = string(parameters(EP)[i])
            pval = real(record_fibre(O)[2])[i]
            summary *= p*raw"""&="""*string(pval)*raw""" \\ """
        end
        summary *= raw"""\end{align*}"""*"\n"
        summary *= raw"""Figure \ref{fig:near_record_fibre} shows a 2-dimensional slice of the parameter space near this fibre."""
        summary *= raw"""\begin{figure}[!htpb]
\centering
\includegraphics[scale=0.5]{OutputFiles/near_record_fibre.png}
\caption{A 2-dimensional slice of the parameter space near the fibre with maximum number of real solutions.}
\label{fig:near_record_fibre}
\end{figure}"""

    (VSD,P) = visualize(EP; near = record_parameters(O), strategy = :quadtree, resolution = 2000)
    save(P, "near_record_fibre.png")
    else
        summary *= raw"""The maximum number of real solutions found in the fibres is $"""*string(maxO)*raw""", but the maximum number of real solutions found in the sample was $"""*string(maxEX)*raw""". """
        summary *= raw"""The fibre for this was not saved."""
    end
    summary *= raw"""The maximum number of real solutions found in the fibres is $"""*string(max_real)*raw"""$. """
    missing_reals = setdiff(collect(possible_reals), EX)
    if in(maxO,missing_reals) == true
        println("!")
        missing_reals = setdiff(missing_reals, maxO)
    end
    if length(missing_reals) == 0
        summary *= raw"""All possible counts of real solutions were found."""
    else
        summary *= raw"""The following counts of real solutions were not found: """
        summary *= raw"""\begin{center}"""*"\n"
        summary *= join(missing_reals, ", ")
        summary *= raw""". """
        summary *= raw"""\end{center}"""*"\n"
    end
end