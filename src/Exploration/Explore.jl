export
    explore





@doc raw"""
    explore(E::EnumerativeProblem,L; sampler = nothing, real_parameters=true, n_samples=1000)

 Explores the parameter space of $E$ while collecting data about the values of the functions listed in $L$ on the fibres. 
   Specifically, $L$ is a list of functions of the form `l(E::EnumerativeProblem,PathRes::Result,Param::Vector{ComplexF64})::Any`
   which are applied to the resulting fibres computed by sampling `n_samples`-many parameters via `sampler`.

 By default, after the sampling and solving procedure, those fibres with fewer than degree(E)-many solutions are thrown out.
  If `real_parameters` is set to `true`, additionally those fibres whose number of real solutions is of opposite parity to degree(E) are thrown out. 
 # Examples
 ```jldoctest
 julia> T = TwentySevenLines();


 julia> Data = explore(T, [real_points_in_fibre,positive_points_in_fibre];n_samples=10000,real_parameters=true)

 Populating a base fibre of the enumerative problem

 Solving for 10000 parameters... 100%|████████████████████████████████████████| Time: 0:00:39
   parameters solved:  10000
   paths tracked:      270000

 Total number of fibres computed:10000
 9937/10000 satisfies degree_check
 9937/10000 satisfies real_parity_check
 2-element Vector{Vector{Any}}:
 ["real_points_in_fibre", "positive_points_in_fibre"]
 [Any[7, 7, 3, 7, 3, 3, 3, 3, 7, 3  …  3, 3, 3, 3, 3, 27, 3, 15, 7, 3], Any[0, 0, 0, 1, 0, 1, 0, 1, 2, 0  …  0, 0, 0, 1, 0, 1, 0, 1, 0, 1]]

 julia> tally(Data[2][1])
 Dict{Any, Int64} with 4 entries:
  3  => 7402
  7  => 2112
  15 => 410
  27 => 13

 julia> tally(Data[2][2])
 Dict{Any, Int64} with 6 entries:
  0 => 7094
  1 => 2283
  2 => 452
  3 => 76
  4 => 27
  5 => 5
 ```
 """
function explore(E::EnumerativeProblem,L; sampler = nothing, real_parameters=true, n_samples=1000)
    Labels = [string(l) for l in L]
    F = system(E)
    P = []
    if sampler!=nothing
        P = sampler(n_samples)
    elseif real_parameters==true
        P = real_sampler(n_samples,n_parameters(E))
    else
        P = complex_sampler(n_samples,n_parameters(E))
    end
    BigSolve = solve_over_params(E,P)
    AllData = []
    for l in L
        Data = []
        for B in BigSolve
            push!(Data,l(E,B[1],B[2]))
        end
    push!(AllData,Data)
    end
    return([Labels,AllData])
end
