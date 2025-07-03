export 
    automate!

"""
    automate!(EP::EnumerativeProblem)

    Automatically applies all registered algorithms in `ALGORITHM_DATA` to the given `EnumerativeProblem` `EP`.
    For each algorithm, the output property is computed for `EP`. 
    This is useful for precomputing or populating all knowable properties of an 
    enumerative problem in a single call.


# Arguments
- `EP::EnumerativeProblem`: The enumerative problem to which all algorithms will be applied.

# Example
```julia
T = TwentySevenLines();
automate!(T)
```
"""
function automate!(EP::EnumerativeProblem)
    for AD in keys(ALGORITHM_DATA) #This scrolls through all algorithms and applies each to EP. 
        if AD != user_given && AD != conjunction
        D = ALGORITHM_DATA[AD]
        ep = output_property(D)
        @vprintln("Computing ",ep, " via ", name(AD))
        ep(EP)
        end
    end
end