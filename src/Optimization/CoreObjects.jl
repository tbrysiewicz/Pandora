export
    Score,
    RealScoreSpace,
    MaxScore


struct Score
    ScoreFunction #must accept Tuple{Result, Vector{ComplexF64}} and return anything in T
    gt #must accept two elements of T and return true/false true = gt or eq
end

function CtsReal(RP::Tuple{Result, Vector{Float64}})
    return((nreal(RP[1]),NonRealMinNorm(RP[1])))
end

function NonRealMinNorm(R::Result)
    record = nothing
    for r in R
        if is_real(r)==false
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

RealScoreSpace = Score(CtsReal,gt)

function MaxScore(Sols::Vector{Tuple{Result,Vector{Float64}}}, SC::Score)
    RecordFibre = Sols[1]
    Record = SC.ScoreFunction(RecordFibre)
    for i in 2:length(Sols)
        newscore=SC.ScoreFunction(Sols[i])
        if SC.gt(newscore,Record)
            RecordFibre=Sols[i]
            Record = newscore
        end
    end
    return((Record,RecordFibre))
end
