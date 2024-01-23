# PANDO(RA)
### Parallel, Automated, Numerical, Discovery and Optimization (Research Aid)

This is a prototype repository for a suite of numerical methods for 
automatically and reliably analyzing enumerative problems. 

### Example: Computing the 27 lines on a cubic surface
 ```julia
julia> using Pandora
julia> E = TwentySevenLines()


           X := V(f_1..f_4) ⊆ C^4 x C^20
           |
           |
           | π    ???-to-1
           |
           V
          C^20

An enumerative problem in 4 variable(s) cut out by 4 condition(s) over 20 parameters.


julia> degree(E)
Populating a base fibre of the enumerative problem
Tracking 45 paths... 100%|██████████| Time: 0:00:03
  # paths tracked:                  45
  # non-singular solutions (real):  27 (0)
  # singular endpoints (real):      0 (0)
  # total solutions (real):         27 (0)
27
```

## Computing the Galois group of the 27 lines
```julia
julia> G = galois_group(E)
Sampling 20 elements from the monodromy group of 

           X := V(f_1..f_4) ⊆ C^4 x C^20
           |
           |
           | π   27-to-1
           |
           V
          C^20

An enumerative problem in 4 variable(s) cut out by 4 condition(s) over 20 parameters.

Tracking 27 paths... 100%|██████████| Time: 0:00:01
  # paths tracked:                  27
  # non-singular solutions (real):  27 (0)
  # singular endpoints (real):      0 (0)
  # total solutions (real):         27 (0)
20 out of 20 are plausible permutations
<permutation group with 20 generators>

julia> using Oscar; order(G)
51840

julia> isprimitive(G)
true

```
### Observing the reality constraints on the 27 lines
```julia
Data = explore(E,[real_points_in_fibre]; n_samples=10000);
Solving for 10000 parameters... 100% Time: 0:00:33
  # parameters solved:  10000
  # paths tracked:      270000
Total number of fibres computed:10000
9965/10000 satisfies degree_check
9965/10000 satisfies real_parity_check


julia> tally(Data[2][1])
Dict{Any, Int64} with 4 entries:
  3  => 7488
  7  => 2059
  15 => 413
  27 => 5
```

The main datatypes in Pandora are currently
 Enumerative Problem
 Variety

The main functionality in Pandora currently includes
- Monodromy Group Computation
- Lazy Evaluation of Enumerative Problems
- Automated Exploration Tools
- (Decorated) Algebraic Matroids

Coming Soon....
- Ability to optimize functions on fibres of enumerative problems over their parameter spaces via shotgun hill-climbing, including
  - Number of real solutions
  - Number of positive solutions
  - Any other upper semicontinuous function
