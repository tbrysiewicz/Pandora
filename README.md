
# PANDO(RA)

## Parallel, Automated, Numerical, Discovery and Optimization (Research Aid)

Pandora is a Julia package for studying enumerative problems with numerical and symbolic methods.

![Pandora logo](ReadMeImages/Pandoralogo.png?raw=true "Parallel, Automated, Numerical, Discovery and Optimization (Research Aid)")

## Contact

Pandora.jl is under active development (current README update: April 2026). Comments and feature requests are welcome.

- Taylor Brysiewicz: tbrysiew@uwo.ca

## Installation

From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

## Quick Start Demo

```julia
using Pandora
using Oscar

@var x, y
@var a[1:4], b[1:4]

f1 = a[1]*x + a[2]*y + a[3]*x^2*y + a[4]*x*y^2
f2 = b[1]*x + b[2]*y + b[3]*x^2*y + b[4]*x*y^2

E = EnumerativeProblem([f1, f2], variables=[x, y], parameters=vcat(a, b), torus_only=true)

degree(E)
bkk_bound(E)
bezout_bound(E)

knowledge(E)
last(knowledge(E))

G = monodromy_group(E)
is_primitive(G)
order(G)
gens(G)

is_decomposable(E)
is_lacunary(E)
```

Expected values for this example include `degree(E) == 4`, `bkk_bound(E) == 4`, and `bezout_bound(E) == 9`.

## Sampling and Optimization

```julia
ER = tally.(explore(E, [n_real_solutions, n_positive_solutions], n_samples=1000))
O = maximize_n_real_solutions(E)
C = certify(record_fibre(O), E)

ER
C
```

## Visualization

```julia
(V, P) = visualize(E; near=record_parameters(O), strategy=:quadtree)
save(P, "MyDiscriminant.png")

refine!(V)
refine!(V)
P2 = visualize(V)
save(P2, "MyBetterDiscriminant.png")

SupportVisualization = visualize_support(E)
save(SupportVisualization[1], "support1.png")
save(SupportVisualization[2], "support2.png")

NP = newton_polytopes(E)
mv = volume(sum(NP)) - volume(NP[1]) - volume(NP[2])
```

Example outputs:

- ![Discriminant](ReadMeImages/MyDiscriminant.png?raw=true "n_real_solutions Visualization")
- ![Refined discriminant](ReadMeImages/MyBetterDiscriminant.png?raw=true "Better n_real_solutions Visualization")
- ![Support visualization](ReadMeImages/support1.png?raw=true "Newton Polytope and Support Visualization")

## Automation and Summary

```julia
T = TwentySevenLines()
automate!(T)
summarize(T)
```

Typical output file:

- `OutputFiles/latex_summary.pdf`

Example summary:

![Pandora summary example](ReadMeImages/PandoraSummaryExample.png?raw=true "Pandora Summary")

## Verification Script (run from repository root)

The following script is a reproducible smoke test for all workflows shown above.

```bash
julia --project=. <<'JULIA'
using Pkg
Pkg.instantiate()

using Pandora
using Oscar

@var x, y
@var a[1:4], b[1:4]

f1 = a[1]*x + a[2]*y + a[3]*x^2*y + a[4]*x*y^2
f2 = b[1]*x + b[2]*y + b[3]*x^2*y + b[4]*x*y^2
E = EnumerativeProblem([f1, f2], variables=[x, y], parameters=vcat(a, b), torus_only=true)

@assert degree(E) == 4
@assert bkk_bound(E) == 4
@assert bezout_bound(E) == 9

k = knowledge(E)
@assert !isempty(k)
@assert last(k) !== nothing

G = monodromy_group(E)
@assert order(G) == 8
@assert length(gens(G)) == 2
@assert is_decomposable(E)
@assert is_lacunary(E)

ER = tally.(explore(E, [n_real_solutions, n_positive_solutions], n_samples=200))
@assert length(ER) == 2

O = maximize_n_real_solutions(E)
C = certify(record_fibre(O), E)
@assert C isa CertificationResult

(V, P) = visualize(E; near=record_parameters(O), strategy=:quadtree)
save(P, "MyDiscriminant.png")

refine!(V)
refine!(V)
P2 = visualize(V)
save(P2, "MyBetterDiscriminant.png")

SupportVisualization = visualize_support(E)
save(SupportVisualization[1], "support1.png")
save(SupportVisualization[2], "support2.png")

NP = newton_polytopes(E)
@assert length(NP) == 2
mv = volume(sum(NP)) - volume(NP[1]) - volume(NP[2])
println("Mixed volume: ", mv)

T = TwentySevenLines()
automate!(T)

try
    summarize(T)
catch err
    @warn "summarize(T) failed (often due to missing LaTeX tools)" exception=(err, catch_backtrace())
end

println("README smoke test completed.")
JULIA
```

## Credits

- Developers
  - Taylor Brysiewicz
  - Alexandra Makris
  - Samuel Feldman
  - Noah Vale
  - Deepak Mundayurvalappil Sadanand