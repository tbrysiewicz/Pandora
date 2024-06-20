export
    tally,
    mean,
    mode,
    stat_print



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

function mode(S)
    T = tally(S)
    K = collect(keys(T))
    FreqMax = max(collect(values(T))...)
    Kmax = filter(x->T[x]==FreqMax,K)
    return(Kmax)
end

function stat_print(S)
    str = "N:"
    str = str*string(length(S))*"\n"
    str = str*"Mean:"*string(mean(S))*"\n"
    str = str*"Mode:"*string(mode(S))*"\n"
    str = str*"Max:"*string(max(S...))*"\n"
    str = str*"Min:"*string(min(S...))*"\n"
    println("N:", length(S))
    println("Mean:", mean(S))
    println("Mode:", mode(S))
    println("Max:", max(S...))
    println("Min:", min(S...))
    return([length(S),mean(S),mode(S),max(S...),min(S...)])
end