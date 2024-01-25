export
    Score,
    RealScoreSpace,
    max_score
   # gt


#TODO: Change style so that commands are consistent w/ rest of Pandora
#      (Structures are Capitalized like EnumerativeProblem, functions are 
#       not, and use underscores like galois_group(E::EnumerativeProblem))
#      Recode the optimization procedures anyways to be more clever as discussed. 
#    
#

 
struct Score
    ScoreFunction #must accept Tuple{Result, Vector{ComplexF64}} and return anything in T
    gt #must accept two elements of T and return true/false true = gt or eq
end


function cts_real(RP::Tuple{Result, Vector{Float64}})
    return((nreal(RP[1]),non_real_min_norm(RP[1])))
end


function non_real_min_norm(R::Result)
    record = nothing                #Here, record ans record_fibre are just variable names, not parts of the struct OptimizerData.
    for r in R
        if HomotopyContinuation.is_real(r)==false
            if record!= nothing
                nm = norm(imag(r.solution))
                if nm<record
                    record=nm
                end
            else
                record = norm(imag(r.solution))
            end
        end
    end
    if record==nothing
        return(0.0)
    else
        return(record)
    end
end

@doc raw"""
    gt(a::Tuple{Int64,Float64},b::Tuple{Int64,Float64})
    
    Compares two elements of the type Tuple{Int64,Float64}, and returns true 
    if the first one is greater than or equal to the latter. 
"""
function gt(a::Tuple{Int64,Float64},b::Tuple{Int64,Float64})
    if a[1]>b[1]
        return(true)
    elseif b[1]>a[1]
        return(false)
    else
        if a[2]<b[2]
            return(true)
        elseif a[2]>b[2]
            return(false)
        else
            return(true)
        end
    end
end

RealScoreSpace = Score(cts_real,gt)

function max_score(Sols::Vector{Tuple{Result,Vector{Float64}}}, SC::Score)
    record_fibre = Sols[1]
    record = SC.ScoreFunction(record_fibre)
    for i in 2:length(Sols)
        newscore=SC.ScoreFunction(Sols[i])
        if SC.gt(newscore,record)
            record_fibre=Sols[i]
            record = newscore
        end
    end
    return((record,record_fibre))
end

