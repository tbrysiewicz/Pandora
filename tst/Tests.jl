TSL = TwentySevenLines()
degree(TSL)
base_fibre(TSL)
base_solutions(TSL)
base_parameters(TSL)


solve_over_params(TSL,[randn(Float64,n_parameters(TSL)) for i in 1:100])

galois_group(TSL)
monodromy_group(TSL)

EP = restrict_enumerative_problem(TSL,[randn(Float64,20) for i in 1:3])

optimize_real(TSL)