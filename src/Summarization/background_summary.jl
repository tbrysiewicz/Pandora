
#########BACKGROUND SUMMARY####################

function texbf(s::String)
    return raw"\textbf{"*s*raw"}"
end

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
    if k==1
        summary *= "\n\n The polynomial system given involves "*raw"1 polynomials $F=(f_"*kk*raw")$ in "
    elseif k==2
        summary *= "\n\n The polynomial system given involves "*raw"2 polynomials $F=(f_1,f_"*kk*raw")$ in "
    elseif k==3
        summary *= "\n\n The polynomial system given involves "*raw"3 polynomials $F=(f_1,f_2,f_"*kk*raw")$ in "
    else
        summary *= "\n\n The polynomial system given involves "*kk*raw" polynomials $F=(f_1,\ldots ,f_"*kk*raw")$ in "
    end
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
          summary *=texbf(string(P[i]))*", "
       end
       summary *=texbf(string(last(P)))*raw").\]"
    else
        summary *= texbf(string(P[1]))*","*texbf(string(P[2]))*raw",\ldots ,"*texbf(string(last(P)))*raw").\]"
    end

    summary *= "The equations look as follows."
    summary *= latex_polynomials(system(EP);ff = "f", kwargs...)
    summary *= raw"\noindent The inequations look as follows."
    summary *= latex_polynomials(inequations(EP); ff = "g", inequations = true, bold_coefficients= false, kwargs...)
    
    summary *= "This setup can be summarized in the following diagram describing the parametrized polynomial system "
    summary *= raw"as a (possibly reducible) branched cover $\pi$ from the incidence variety $\mathcal V(F)$ to the parameter space $\mathbb{C}^"*NN*raw"$."
    summary *= tikz_incidence(EP)


    return(summary)
end
function tikz_incidence(EP::EnumerativeProblem)
    n = ambient_dimension(EP)
    k = n_polynomials(EP)
    kiq = n_inequations(EP)
    N = n_parameters(EP)

    nn = raw"{"*string(n)*raw"}"
    kk = raw"{"*string(k)*raw"}"
    kkiq = raw"{"*string(kiq)*raw"}"
    NN = raw"{"*string(N)*raw"}"

    s = raw"""
    \begin{center}
    \begin{tikzpicture}
    \matrix (m) [matrix of math nodes,row sep=3em,column sep=0em,minimum width=2em]
    {
        \mathcal V(F)"""
    if length(inequations(EP))>0
    s = s*raw"""-  \mathcal V(G)"""
    end
    s=s*raw"""=\mathcal V("""
    if k == 1
        s *= raw"f_1"
    elseif k == 2
        s *= raw"f_1,f_2"
    elseif k == 3
        s *= raw"f_1,f_2,f_3"
    else
        s *= raw"f_1,\ldots,f_" * kk * raw""
    end
    s *= raw""")"""
    if kiq>0
        s *= raw""" - \mathcal V("""
        if kiq == 1
            s *= raw"g_1"
        elseif kiq == 2
            s *= raw"g_1,g_2"
        elseif kiq == 3
            s *= raw"g_1,g_2,g_3"
        else
            s *= raw"_1,\ldots,g_" * kkiq 
        end
        s *= raw""")"""
    end
    s *= raw"""=\{(\textbf{x},\textbf{p}) \mid F(\textbf{x};\textbf{p})=\textbf{0}"""
    if kiq>0
        s *= raw""", G(\textbf{x}) \neq 0"""
    end
    s *= raw"""\} & \subset \mathbb{C}_{\textbf{x}}^""" * nn * raw""" \times \mathbb{C}_{\textbf{p}}^""" * NN * raw""" \\
        \mathbb{C}_{\p}^""" * NN * raw""" &  \\};
    \path[-stealth]
        (m-1-1) edge node [left] {$\pi$} (m-2-1);
    \end{tikzpicture}
    \end{center}
    """
    return s
end

function monomial_string(BasePowerPairs)
    ms = raw""
    for (v,p) in BasePowerPairs
        if p != 0
            if p ==1
                ms *= string(v)
            else
                ms *= string(v)*raw"^{ "*string(p)*raw"}"
            end
        end
    end
    return(ms)
end

function terms_as_text_list(f::Expression,V::Vector{Variable},P::Vector{Variable}; bold_coefficients = true)
    F = System([f],variables=V,parameters=P)
    SC = support_coefficients(F)
    M = SC[1][1]
    C = SC[2][1]
    (n,n_terms)=size(M)
    if bold_coefficients
        myterms = reverse([raw"{\textbf{"*string(C[i])*raw"}}"*raw" \cdot "*monomial_string([(V[j],M[j,i]) for j in 1:n]) for i in 1:n_terms])
        return(myterms)
    else
        C = map(c->c==1 ? "" : string(c)*raw"\cdot",C)
        myterms = reverse([raw"{{"*string(C[i])*raw"}}"*monomial_string([(V[j],M[j,i]) for j in 1:n]) for i in 1:n_terms])
        return(myterms)
    end
end

function latex_polynomials(FF::System;ff = "f", inequations = false, kwargs...)
    V = variables(FF)
    P = parameters(FF)
    F = expressions(FF)
    s = raw"$$ \begin{array}{crc}"*"\n"
    for i in 1:length(F)
        f = F[i]
        t = terms_as_text_list(f,V,P; kwargs...)
        s *= ff*raw"_"*string(i)*raw": &"
        if length(t)<10
            polystring = t[1]
            for j in 2:length(t)
                polystring *= raw" + "*t[j]
            end
            s *= replace(polystring,"*"=>"")
        else
            f_start = t[1]
            f_start *= "+"*t[2]*"+"*t[3]*"+"
            f_end = t[end-1]*"+"*t[end]
            f_replaceable = f_start*raw"\cdots "*f_end

            s *= replace(string(f_replaceable),"*"=>"")
        end
        if inequations == true
            s *= raw"&\neq 0"
        else
         s *= raw"&= 0"
        end
        if i<length(F)
            s*=raw"""\\\\"""*"\n"
        end
    end
    s *= "\n"*raw"\end{array}$$"*"\n\n"
    return(s)
        
end
