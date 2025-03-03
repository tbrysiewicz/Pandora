export 
    automate!

function automate!(EP::EnumerativeProblem)
    for AD in keys(ALGORITHM_DATA)
        D = ALGORITHM_DATA[AD]
        ep = output_property(D)
        println("Computing ",ep)
        ep(EP)
    end
end