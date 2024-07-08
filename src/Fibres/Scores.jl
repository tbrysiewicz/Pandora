export
    Score,
    TrivialScore,
    score_function


struct Score
    ScoreFunction::Function #must accept Tuple{Result, Vector{ComplexF64}} and return anything in T
    gt::Function #must accept two elements a,b of T and return true/false, true means a>=b
end

function score_function(Score)
    Score.ScoreFunction
end

const TrivialScore() = Score(RP -> 0,<)