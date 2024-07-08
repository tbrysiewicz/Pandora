import Oscar.is_decomposable
import Oscar.galois_group

export
    galois_group,
    is_decomposable,
    monodromy_group

@doc raw"""
    galois_group(E::EnumerativeProblem; nloops = 20)

 Computes the Galois group, or monodromy group, of an enumerative problem. 

 `nloops` indicates how many monodromy loops are taken to generate the Galois group


 # Examples
 ```jldoctest

 julia> T = TwentySevenLines();

 julia> G = galois_group(T)
 Sampling 20 elements from the monodromy group of 

           X := V(f_1..f_4) ⊆ C^4 x C^20
           |
           |
           | π    ???-to-1
           |
           V
          C^20

 An enumerative problem in 4 variable(s) cut out by 4 condition(s) over 20 parameters.

 Populating a base fibre of the enumerative problem
 19 out of 20 are plausible permutations
 <permutation group with 19 generators>

 julia> Oscar.order(G)
 51840
 ```
 """
function galois_group(E::EnumerativeProblem; nloops = 20, radius=1)
    S = sample_monodromy_elements(E,nloops; radius=radius)
    return(permutation_group(degree(E),S))
end

@doc raw"""
    monodromy_group(E::EnumerativeProblem; nloops = 20)

 Computes the Galois group, or monodromy group, of an enumerative problem. 

 `nloops` indicates how many monodromy loops are taken to generate the Galois group


 # Examples
 ```jldoctest

 julia> T = TwentySevenLines();

 julia> G = monodromy_group(T)
 Sampling 20 elements from the monodromy group of 

           X := V(f_1..f_4) ⊆ C^4 x C^20
           |
           |
           | π    ???-to-1
           |
           V
          C^20

 An enumerative problem in 4 variable(s) cut out by 4 condition(s) over 20 parameters.

 Populating a base fibre of the enumerative problem
 20 out of 20 are plausible permutations
 <permutation group with 20 generators>

 julia> Oscar.order(G)
 51840
 ```
 """
function monodromy_group(E::EnumerativeProblem; nloops = 20, radius = 1)
    S = sample_monodromy_elements(E,nloops; radius = radius)
    return(permutation_group(degree(E),S))
end

@doc raw"""
    is_decomposable(E::EnumerativeProblem)

 Establishes if an enumerative problem is decomposable by checking if its monodromy group is imprimitive. 
 """
function is_decomposable(E::EnumerativeProblem)
    G = galois_group(E)
    return Oscar.isprimitive(G)==false
end
