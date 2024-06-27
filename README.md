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

### Finding all real solutions
```julia
julia> OD = optimize_enumerative(TSL,RealScoreSpace,3;bucket_size=500)
Nparameters:20
Total number of fibres computed:1
1/1 satisfies degree_check
1/1 satisfies real_parity_check
Step: 1
Old Radius: 1.0
Radius: 1.0
Total number of fibres computed:503
501/503 satisfies degree_check
501/503 satisfies real_parity_check
Progress:   0.22613184080819235
Record:(15, 0.24857517194116546)
Taboo Score:0.6866267465069861
Closeness of the real solutions: 1.0230028380217617
Stuck Score:0
-----------------------------------------------------------
Step: 2
Old Radius: 1.0
Radius: 1.0
Total number of fibres computed:503
501/503 satisfies degree_check
501/503 satisfies real_parity_check
Progress:   0.046864661191914056
Record:(27, 0.0)
Taboo Score:0.8842315369261478
Closeness of the real solutions: 1.141065805355264
Stuck Score:1
-----------------------------------------------------------
Step: 3
Old Radius: 1.0
Radius: 0.1
Total number of fibres computed:503
499/503 satisfies degree_check
499/503 satisfies real_parity_check
Progress:   -0.0
Record:(27, 0.0)
Taboo Score:0.3346693386773547
Closeness of the real solutions: 0.5851193752114505
Stuck Score:2
-----------------------------------------------------------
OptimizerData((Result with 27 solutions
========================
• 27 paths tracked
• 27 non-singular solutions (27 real)
• random_seed: 0x24e8119b
, [-3.0276536401879435, -3.1698278462259983, 0.3285065039454872, 0.5608183584294572, -1.124399411482888, 1.1403466968409988, -0.2233627471248292, 2.5923566620936898, 1.2415820889520897, 1.2945705189804517, 4.277148441687768, -1.2470039966787618, 0.041635894443212695, 3.167801437241401, 0.18925573552102234, 1.642076066998358, -1.3380539376892606, -2.353872405568242, -0.6912648497026905, -1.4426115948474645]), (27, 0.0), (0.3346693386773547, 0.5851193752114505), 2, (Result with 27 solutions
========================
• 27 paths tracked
• 27 non-singular solutions (27 real)
• random_seed: 0xdfd07b4f
, [-2.8241777434824575, -3.027631906368397, 0.25412431237383715, 0.5041502632250112, -1.0105245167569203, 0.967884760082899, -0.16552900658040492, 2.538105981831665, 1.1952005920777247, 1.16431579367716, 4.021897851249387, -1.2477984746846822, -0.11926023635077332, 3.1071980299959434, 0.18124954974875376, 1.6810077739136386, -1.310038511005228, -2.216173074895907, -0.5276682347138648, -1.435545410711053]), 0.1, [0.8409852901473666, 0.8650682697838232, -0.8399031126921687, -0.22741511065499778, 0.04700318358050594, 0.016219619978378924, -1.0798252035357843, -0.7050019400980675, 1.0132929125316694, -0.941780523050524, 0.6030474171265916, -0.8615554917243122, 0.28279089646800026, 0.9876301010470466, -0.25085170131749, -0.6395896913655517, 0.5102548566968205, 1.8742945404959153, -1.3699761880782597, 0.11763816145977886], Pandora.Strategies(false, true, false, false, true), 1)

```

### Visualizing a discriminant slice
```


julia> P = OD.RecordFibre[2]
20-element Vector{Float64}:
 -3.0276536401879435
 -3.1698278462259983
  0.3285065039454872
  0.5608183584294572
 -1.124399411482888
  1.1403466968409988
 -0.2233627471248292
  2.5923566620936898
  1.2415820889520897
  1.2945705189804517
  4.277148441687768
 -1.2470039966787618
  0.041635894443212695
  3.167801437241401
  0.18925573552102234
  1.642076066998358
 -1.3380539376892606
 -2.353872405568242
 -0.6912648497026905
 -1.4426115948474645
 
 
julia> EP = restrict_enumerative_problem(TSL,[P+randn(Float64,20),P+randn(Float64,20),P]);

julia> visualize_discriminant(EP, "Delaunay";resolution=5000,total_resolution=20000,depth=4,scatter=true,xlims=[-2,2],ylims=[-2,2])

```


![Alt text](Discriminant.png?raw=true "Discriminant Visualization")

The main datatypes in Pandora are currently
 Enumerative Problem
 Variety

The main functionality in Pandora currently includes
- Monodromy Group Computation
- Lazy Evaluation of Enumerative Problems
- Automated Exploration Tools
- Ability to optimize functions on fibres of enumerative problems over their parameter spaces via shotgun hill-climbing, including
  - Number of real solutions
  - Number of positive solutions
  - Any other given score function
