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
27

julia> degree(E;Method = "Monodromy",Retry = true)
Populating a base fibre of the enumerative problem
Using monodromy
27
```

## Computing the Galois group of the 27 lines
```julia

julia> G = galois_group(E; nloops = 100, radius=10)
Sampling 100 elements from the monodromy group of 

           X := V(f_1..f_4) ⊆ C^4 x C^20
           |
           |
           | π   27-to-1
           |
           V
          C^20

An enumerative problem in 4 variable(s) cut out by 4 condition(s) over 20 parameters.

99 out of 100 are plausible permutations
Permutation group of degree 27

julia> using Oscar; order(G)
51840

julia> isprimitive(G)
true

```
### Observing the reality constraints on the 27 lines
```julia
julia> Data = explore(E,[real_points_in_fibre];n_samples=10000);
Solving for 10000 parameters... 100%|████████████████████████████████████████| Time: 0:00:14
  # parameters solved:  10000
  # paths tracked:      270000
Total number of fibres computed:10000
9957/10000 satisfies degree_check
Retry:false
LENGTH OF SOLS:9957


julia> tally(Data[3][1])
Dict{Any, Int64} with 4 entries:
  3  => 7488
  7  => 2059
  15 => 413
  27 => 5
```

### Finding all real solutions
```julia
julia> OD = optimize_real(E)
i'i.i'i.!
---------------------------------------------------------------------------
Optimizer for an enumerative problem with 27 many solutions.
  Current barrier-weighted objective function: (27.0, -0.04042909428237266)
   Record barrier-weighted objective function: (27.0, -0.04042909428237266)
  Goal: Reached
---------------------------------------------------------------------------



```

### Visualizing a discriminant slice
```
julia> P = OD.record_fibre[2]
20-element Vector{Float64}:
 -2.804451318458232
  0.7940808007947211
  7.361465112381069
 -4.611297477231954
 -0.27901826818369635
 -0.7923994094742532
 -6.573816057838942
  ⋮
 -0.7657799067649271
  1.3448489447510072
  1.9205988028997956
  5.550237971895067
 -1.0080509935587196
  1.0963234606809722
 
julia> EP = restrict_enumerative_problem(E,[P+randn(Float64,20),P+randn(Float64,20),P]);


julia> (ImageData,Image) = visualize_parameter_space(EP);

```


![Alt text](Discriminant.png?raw=true "Discriminant Visualization")

The main datatypes in Pandora are currently
- Enumerative Problem
- Variety

The main functionality in Pandora currently includes
- Monodromy Group Computation
- Lazy Evaluation of Enumerative Problems
- Automated Exploration Tools
- Ability to optimize functions on fibres of enumerative problems over their parameter spaces via shotgun hill-climbing, including
  - Number of real solutions
  - Number of positive solutions
  - Any other given score function
