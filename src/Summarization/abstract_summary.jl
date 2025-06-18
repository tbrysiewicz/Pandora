export 
    abstract_summary
#########ABSTRACT SUMMARY####################

function abstract_summary(EP::EnumerativeProblem;kwargs...)
    abstract = raw"""This \LaTeX document, created automatically by Pandora.jl, describes several properties of an 
                     enumerative problem."""
    latex_abstract = raw"""\begin{abstract} 
    """*abstract*raw"""\end{abstract}"""
    return(latex_abstract)
end
