
#########BACKGROUND SUMMARY####################

function background_summary(EP::EnumerativeProblem;kwargs...)
    n = ambient_dimension(EP)
    k = n_polynomials(EP)
    N = n_parameters(EP)
    V = variables(system(EP))
    P = parameters(system(EP))

    nn = raw"{"*string(n)*raw"}"
    kk = raw"{"*string(k)*raw"}"
    NN = raw"{"*string(N)*raw"}"

    summary = raw"""\section{Background}"""
    summary *= "\n\n The polynomial system given involves "*kk*raw" polynomials $F=(f_1,\ldots ,f_"*kk*raw")$ in "
    summary *= nn*raw" variables \[ \x = ("
    if n<10
       for i in 1:n-1
          summary *=string(V[i])*", "
       end
       summary *=string(last(V))*raw") \]"*" and "
    else
        summary *= string(V[1])*","*string(V[2])*raw",\ldots ,"*string(last(V))*raw")\] and "
    end
    summary *= NN*raw" parameters \[\p=("
    if N<10
       for i in 1:N-1
          summary *=string(P[i])*", "
       end
       summary *=string(last(P))*raw").\]"
    else
        summary *= string(P[1])*","*string(P[2])*raw",\ldots ,"*string(last(P))*raw").\]"
    end

    summary *= "The polynomials look as follows."
    summary *= latex_polynomials(EP;kwargs...)
    
    summary *= "This setup can be summarized in the following diagram describing the parametrized polynomial system "
    summary *= raw"as a branched cover $\pi$ from the incidence variety $\mathcal V(F)$ to the parameter space $\mathbb{C}^"*NN*raw"$."
    summary *= tikz_incidence(EP)


    return(summary)
end

function tikz_incidence(EP::EnumerativeProblem)
    n = ambient_dimension(EP)
    k = n_polynomials(EP)
    N = n_parameters(EP)

    nn = raw"{"*string(n)*raw"}"
    kk = raw"{"*string(k)*raw"}"
    NN = raw"{"*string(N)*raw"}"

s = raw"""
    \begin{center}
    \begin{tikzpicture}
    \matrix (m) [matrix of math nodes,row sep=3em,column sep=0em,minimum width=2em]
    {
        \mathcal V(F)=\mathcal V("""*raw"""f_1,\ldots,f_"""*kk*raw""")=\{(\textbf{x},\textbf{p}) \mid F(\textbf{x};\textbf{p})=\textbf{0}\} & \subset \mathbb{C}_{\textbf{x}}^"""*nn*raw""" \times \mathbb{C}_{\textbf{p}}^"""*NN*raw""" \\
        \mathbb{C}_{\p}^"""*NN*raw""" &  \\};
    \path[-stealth]
        (m-1-1) edge node [left] {$\pi$} (m-2-1);
    \end{tikzpicture}
    \end{center}
    """
    return(s)
end

function terms_as_text_list(f::Expression,V::Vector{Variable},P::Vector{Variable})
    F = System([f],variables=V,parameters=P)
    SC = support_coefficients(F)
    M = SC[1][1]
    C = SC[2][1]
    (n,n_terms)=size(M)

    myterms = reverse([raw"{\textbf{"*string(C[i])*raw"}}"*raw" \cdot "*string(prod([V[j]^M[j,i] for j in 1:n])) for i in 1:n_terms])
    return(myterms)
end

function latex_polynomials(EP::EnumerativeProblem;kwargs...)
    V = variables(system(EP))
    P = parameters(system(EP))
    F = expressions(system(EP))
    s = raw" \begin{align*}"*"\n"
    for i in 1:length(F)
        f = F[i]
        t = terms_as_text_list(f,V,P)
        if length(t)<10
            s *= replace(string(f),"*"=>"")
        else
            f_start = t[1]
            f_start *= "+"*t[2]*"+"*t[3]*"+"
            f_end = t[end-1]*"+"*t[end]
            f_replaceable = f_start*raw"\cdots "*f_end

            s *= replace(string(f_replaceable),"*"=>"")
        end
         s *= raw"&= 0"
        if i<length(F)
            s*=raw"""\\\\"""*"\n"
        end
    end
    s *= "\n"*raw"\end{align*}"*"\n\n"
    return(s)
        
end
