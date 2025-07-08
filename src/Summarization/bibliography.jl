export 
    bibliography_summary

function bibliography_summary(EP::EnumerativeProblem; kwargs...)
    summary = raw"\section{Bibliography}"
    
    # Add the bibliography section
    summary *= raw"""TBD"""
end