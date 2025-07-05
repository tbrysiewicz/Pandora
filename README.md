
# PANDO(RA)
### Parallel, Automated, Numerical, Discovery and Optimization (Research Aid)

Pandora is code for studying enumerative problems using numerical and symbolic methods. 



![Alt text](Pandoralogo.png?raw=true "Parallel, Automated, Numerical, Discovery and Optimization (Research Aid)")




## Quick Start Demo

Below is a typical workflow, with code and expected output. 

---

### 1. **Define Variables and Parameters**

```julia
@var x, y
@var a[1:4], b[1:4]
```

### 2. **Define Your System**

```julia
f1 = a[1]*x + a[2]*y + a[3]*x^2*y + a[4]*x*y^2
f2 = b[1]*x + b[2]*y + b[3]*x^2*y + b[4]*x*y^2
E = EnumerativeProblem([f1,f2],variables = [x,y], parameters = vcat(a,b), torus_only=true)


           X := V(f_1..f_2) ⊆ C^2 x C^8
           |
           |
           | π   4-to-1
           |
           V
          C^8

An enumerative problem in 2 variable(s) cut out by 2 condition(s) over 8 parameter(s).

julia> degree(E)
4

julia> bkk_bound(E)
4

julia> bezout_bound(E)
9
```
### 3. **Information is collected in knowledge field of the enumerative problem**

```julia
julia> knowledge(E)
1) [system] as computed by (User given) applied to (nothing).
2) [inequations] as computed by (User given) applied to (nothing).
3) [base_fibre] as computed by (polyhedral_homotopy) applied to [system, inequations]
4) [degree] as computed by (n_solutions) applied to [base_fibre]
5) [bkk_bound] as computed by (bkk_bound) applied to [system]
6) [degree_sequence] as computed by (degree_sequence) applied to [system]
7) [bezout_bound] as computed by (bezout_bound) applied to [degree_sequence]


julia> last(knowledge(E))
└── [bezout_bound] as computed by (bezout_bound)
    └── [degree_sequence]
        └── [system]
```


### 4. **Explore Monodromy Groups as Oscar Objects**

```julia
G = monodromy_group(E)
Permutation group of degree 4

is_primitive(G)
false

order(G)
8

gens(G)
2-element Vector{Oscar.PermGroupElem}:
 (3,4)
 (1,3)(2,4)

is_decomposable(E)
true

is_lacunary(E)
true
```

---

### 5. **Sample From the Parameter Space and Optimize**

```julia
ER = tally.(explore(E, [n_real_solutions, n_positive_solutions]), n_samples = 1000)
Number of valid fibres:1000
2-element Vector{Dict{Any, Int64}}:
 Dict(0 => 426, 4 => 249, 2 => 325)
 Dict(0 => 624, 2 => 23, 1 => 353)

O = maximize_n_real_solutions(E)
C = certify(record_fibre(O),E)
CertificationResult
===================
• 4 solution candidates given
• 4 certified solution intervals (4 real, 0 complex)
• 4 distinct certified solution intervals (4 real, 0 complex)

```

---

### 6. **Visualize and Save Discriminant**

```julia
(V,P) = visualize(E; near = record_parameters(O), strategy = :quadtree)
EP consists of more than two parameters. Visualizing a random 2-plane in the parameter space.
Resolution used:711
Resolution used:1498
save(P, "MyDiscriminant.png")
```

![Alt text](OutputFiles/MyDiscriminant.png?raw=true "n_real_solutions Visualization")


### 7. **Refine and Save Improved Visualization**

```julia
refine!(V);
refine!(V);
P = visualize(V)
save(P, "MyBetterDiscriminant.png")
```

![Alt text](OutputFiles/MyBetterDiscriminant.png?raw=true "Better n_real_solutions Visualization")


### 8. **Visualize Supports and Newton Polytopes**

```julia
SupportVisualization = visualize_support(E);
save(SupportVisualization[1], "support1.png");
save(SupportVisualization[2], "support2.png");
```

![Alt text](OutputFiles/support1.png?raw=true "Newton Polytope and Support Visualization")
![Alt text](OutputFiles/support2.png?raw=true "Newton Polytope and Support Visualization")


```julia

NP = newton_polytopes(E)
2-element Vector{Oscar.Polyhedron}:
 Polytope in ambient dimension 2
 Polytope in ambient dimension 2

mv = volume(sum(NP)) - volume(NP[1]) - volume(NP[2])   # 4
```

### 9. **Automate Knowledge and Summarize**

```julia

julia> T = TwentySevenLines()


           X := V(f_1..f_4) ⊆ C^4 x C^20
           |
           |
           | π   27-to-1
           |
           V
          C^20

An enumerative problem in 4 variable(s) cut out by 4 condition(s) over 20 parameter(s).
The following information is known about this problem:
-system
-inequations
-base_fibre
-degree


julia> automate!(T)
Computing newton_polytopes via newton_polytopes
Pandora.jl is automatically finding an algorithm to compute newton_polytopes. To specify an algorithm, call again with algorithm=>[nameofalgorithm]
There is a total of 1 algorithm(s) in Pandora.jl which compute(s) newton_polytopes:
      1) [USING] newton_polytopes
Computing bezout_bound via bezout_bound
Pandora.jl is automatically finding an algorithm to compute bezout_bound. To specify an algorithm, call again with algorithm=>[nameofalgorithm]
There is a total of 1 algorithm(s) in Pandora.jl which compute(s) bezout_bound:
      1) [USING] bezout_bound
Pandora.jl is automatically finding an algorithm to compute degree_sequence. To specify an algorithm, call again with algorithm=>[nameofalgorithm]
There is a total of 1 algorithm(s) in Pandora.jl which compute(s) degree_sequence:
      1) [USING] degree_sequence
Computing base_fibre via polyhedral_homotopy
Computing degree via n_solutions
Computing degree_sequence via degree_sequence
Computing monodromy group via group_generated_by_monodromy_loops
Pandora.jl is automatically finding an algorithm to compute monodromy group. To specify an algorithm, call again with algorithm=>[nameofalgorithm]
There is a total of 1 algorithm(s) in Pandora.jl which compute(s) monodromy group:
      1) [USING] group_generated_by_monodromy_loops
Pandora.jl is automatically finding an algorithm to compute monodromy_sample. To specify an algorithm, call again with algorithm=>[nameofalgorithm]
There is a total of 1 algorithm(s) in Pandora.jl which compute(s) monodromy_sample:
      1) [USING] Sample of [n_monodromy_loops] random loops of scaling [monodromy_loop_scaling]
# Loops computed:            50
# Valid permutations:        50
# Unique valid permutations: 50
Computing bkk_bound via bkk_bound
Pandora.jl is automatically finding an algorithm to compute bkk_bound. To specify an algorithm, call again with algorithm=>[nameofalgorithm]
There is a total of 1 algorithm(s) in Pandora.jl which compute(s) bkk_bound:
      1) [USING] bkk_bound
Computing affine_bkk_bound via affine_bkk_bound
Pandora.jl is automatically finding an algorithm to compute affine_bkk_bound. To specify an algorithm, call again with algorithm=>[nameofalgorithm]
There is a total of 1 algorithm(s) in Pandora.jl which compute(s) affine_bkk_bound:
      1) [USING] affine_bkk_bound
Computing support via support
Pandora.jl is automatically finding an algorithm to compute support. To specify an algorithm, call again with algorithm=>[nameofalgorithm]
There is a total of 1 algorithm(s) in Pandora.jl which compute(s) support:
      1) [USING] support
Computing fibre_datum via Fibre Datum
Pandora.jl is automatically finding an algorithm to compute fibre_datum. To specify an algorithm, call again with algorithm=>[nameofalgorithm]
There is a total of 1 algorithm(s) in Pandora.jl which compute(s) fibre_datum:
      1) [USING] Fibre Datum
Computing monodromy_sample via Sample of [n_monodromy_loops] random loops of scaling [monodromy_loop_scaling]

julia> summarize(T)
This is pdfTeX, Version ......
..............................
```
**Output:**  
- Summary written to `OutputFiles/latex_summary.pdf`


### 10. An example of the kind of output on the enumerative problem E above


![Alt text](PandoraSummaryExample.png?raw=true "Pandora Summary")

