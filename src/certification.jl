export 
    certify_n_real,
    n_solutions_certified

function certify_n_real(EP::EnumerativeProblem, F::Fibre)
    C = certify(system(EP),F[1]; target_parameters = F[2])
    n_certified = count(x->x.certified, C.certificates)
    n_real_certified = count(x->x.real,C.certificates)
    if n_certified == degree(EP)
        return(n_real_certified)
    else
        return(nothing)
    end
end



function n_solutions_certified(EP::EnumerativeProblem, F::Fibre)
    C = certify(system(EP),F[1]; target_parameters = F[2])
    return(count(x->x.certified, C.certificates))
end