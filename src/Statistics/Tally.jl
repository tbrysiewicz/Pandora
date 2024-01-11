export
    tally,
    mean



@doc raw"""
    tally(S::Vector{Any})

 Produces a dictionary whose keys are the unique elements of `S` and whose value at a key `s` is the number of instances of `s` in `S`.
 """
function tally(S::Vector{Any})
    D=Dict{Any,Int}()
    for s in S
        D[s]=get(D,s,0)+1
    end
    return D
end

function mean(S)
    return(sum(S)/length(S))
end
