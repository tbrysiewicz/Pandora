export 
    certify_n_real,
    n_solutions_certified,
    certify




function is_certified(FD::FibreDatum)
    count(x->x.certified, certificates(FD))==length(FD.solutions)
end

has_certification_result(FD::FibreDatum) = !isnothing(FD.certificates)
certification_result(FD::FibreDatum) = FD.certificates
certificates(FD::FibreDatum) = FD.certificates.certificates
system(FD::FibreDatum) = FD.F

function certify!(FD::FibreDatum) 
    C = certify(F, S; target_parameters = P)
    FD.certificates = C
end

function certify(FD::FibreDatum)
    # Ensure that the system and solutions are set
    if isnothing(FD.F) 
        @error "FibreDatum must have a system before certification. Call `certify(FD::FibreDatum,F::System)` or `certify(FD::FibreDatum, EP::EnumerativeProblem)`"
    end
    certify(system(FD), solutions(FD); target_parameters = parameters(FD))
end

certify(FD::FibreDatum,F::System) = certify(F, solutions(FD); target_parameters = parameters(FD))
certify(FD::FibreDatum, EP::EnumerativeProblem) = certify(system(EP), solutions(FD); target_parameters = parameters(FD))


function certify!(FD::FibreDatum,F::System) 
    FD.certificates = certify(F, solutions(FD); target_parameters = parameters(FD))
    FD.F = F
end
function certify!(FD::FibreDatum, EP::EnumerativeProblem)
    FD.certificates = certify(system(EP), solutions(FD); target_parameters = parameters(FD))
    FD.F = system(EP)
end


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