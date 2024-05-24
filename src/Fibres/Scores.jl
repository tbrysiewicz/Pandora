export
    Score,


struct Score
    ScoreFunction #must accept Tuple{Result, Vector{ComplexF64}} and return anything in T
    gt #must accept two elements of T and return true/false true = gt or eq
end
 