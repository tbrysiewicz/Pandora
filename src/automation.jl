export 
    automate!

function automate!(EP::EnumerativeProblem)
    for AD in keys(ALGORITHM_DATA) #This scrolls through all algorithms and applies each to EP. 
        D = ALGORITHM_DATA[AD]
        ep = output_property(D)
        if ep!= ANY
            println("Computing ",ep)
            ep(EP)
        end
    end
end