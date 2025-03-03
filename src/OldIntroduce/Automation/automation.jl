export 
    automate!

function automate!(EP::EnumerativeProblem)
    for EA in MAIN_ALGORITHMS
        ep = output_property(EA)
        println("Computing ",ep)
        ep(EP)
    end
end