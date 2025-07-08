export 
    certify_n_real,
    n_solutions_certified,
    n_real_solutions_certified,
    certify,
    certify!,
    certificates



has_certification_result(FD::FibreDatum) = !isnothing(FD.certificates)

function is_certified(FD::FibreDatum)
    has_certification_result(FD) ? count(x->x.certified, certificates(FD))==length(FD.solutions) : false
end

certification_result(FD::FibreDatum) = FD.certificates
certificates(FD::FibreDatum) = FD.certificates.certificates
system(FD::FibreDatum) = FD.F
has_system(F::Fibre) = F isa FibreDatum ? system(F)!=nothing : false


function certify(fibre::Fibre, F::System)::CertificationResult
    certify(F, fibre[1]; target_parameters = fibre[2])
end

function certify(fibre::Fibre, EP::EnumerativeProblem)::CertificationResult
    F = system(EP)
    certify(F, fibre[1]; target_parameters = fibre[2])
end


function certify(F::Fibre) ::CertificationResult
    if has_system(F)
        return(certify(F,system(F)))
    else
        @error "Fibre must have a reference to a system before certification. Call `certify(F::Fibre,F::System)` or `certify(F::Fibre, EP::EnumerativeProblem)`"
    end
end

function certify!(FD::FibreDatum, F::System) ::CertificationResult
    C = certify(FD, F)
    FD.certificates = C
    FD.F = F
    return(C)
end

function certify!(FD::FibreDatum)::CertificationResult
    if has_system(FD)
        return(certify!(FD, system(FD)))
    else
        error("FibreDatum must have a reference to a system before certification. Call `certify!(FD::FibreDatum,F::System)` or `certify!(FD::FibreDatum, EP::EnumerativeProblem)`")
    end
end

function certify!(FD::FibreDatum, EP::EnumerativeProblem)  ::CertificationResult
    return(certify!(FD, system(EP)))
end

function n_solutions_certified(FD::FibreDatum)
    if is_certified(FD)
        return(count(x->x.certified, certificates(FD)))
    else
        return(nothing)
    end
end
function n_real_solutions_certified(FD::FibreDatum)
    if is_certified(FD)
        return(count(x->x.certified && x.real, certificates(FD)))
    else
        return(nothing)
    end
end

#=

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

=#