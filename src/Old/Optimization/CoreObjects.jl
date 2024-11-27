export
    RealScoreSpace,
    RealScoreSpaceTotal,
    RealScoreSpaceNoImag,
    RealScoreMin,
    max_score,
    real_min_dist
   # gt


#TODO: Change style so that commands are consistent w/ rest of Pandora
#      (Structures are Capitalized like EnumerativeProblem, functions are 
#       not, and use underscores like galois_group(E::EnumerativeProblem))
#      Recode the optimization procedures anyways to be more clever as discussed. 
#    
#

function cts_real(RP::Tuple{Result, Vector{Float64}})
    return((nreal(RP[1]),non_real_min_norm(RP[1])))
end

function cts_real_total(RP::Tuple{Result, Vector{Float64}})
    return((nreal(RP[1]),non_real_total_norm(RP[1])))
end

function cts_real_no_imag(RP::Tuple{Result, Vector{Float64}})
    return((nreal(RP[1])-5*n_pure_imag(RP[1]),non_real_min_norm(RP[1])+100000*pure_imag_closeness(RP[1])))
end

function cts_real_min(RP::Tuple{Result,Vector{Float64}})   ##for finding the fewest number of real solutions
    return((-nreal(RP[1]),real_min_dist(RP[1])))
end

function real_min_dist(R::Result)
S = HomotopyContinuation.real_solutions(R)
if length(S)==0
  return(0.0)
end
record=Inf
n = length(S)
for i in 1:n
  for j in 1:n
    if i!=j
       s1=S[i]
       s2=S[j]
       nm = norm(s1-s2)
       if nm<record
          record=nm
       end
    end
   end
end
return(record)
end



function pure_imag_closeness(R::Result)
    PI = filter(x->norm(real(x))<0.0000000001,solutions(R))
    di=-1
    for i in 1:length(PI)
        for j in 1:length(PI)
            if i!=j
                n = norm(PI[i]-PI[j])
                if n<di || di==-1
                    di=n
                end
            end
        end
    end
    if di==-1
        di=0
    end
    return(di)
end


function  n_pure_imag(R::Result)
    length(filter(x->norm(real(x))<0.0000000001,solutions(R)))
end


function non_real_total_norm(R::Result)
    total = 0.0                #Here, record ans record_fibre are just variable names, not parts of the struct OptimizerData.
    for r in R
        if HomotopyContinuation.is_real(r)==false
            total=total+norm(imag(r.solution))
        end
    end
    return(total)
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
RealScoreSpaceTotal = Score(cts_real_total,gt)
RealScoreSpaceNoImag = Score(cts_real_no_imag,gt)
RealScoreMin=Score(cts_real_min,gt)


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

