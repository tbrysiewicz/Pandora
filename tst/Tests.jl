TSL = TwentySevenLines()
degree(TSL)
base_fibre(TSL)
base_solutions(TSL)
base_parameters(TSL)

galois_group(TSL)
monodromy_group(TSL)

EP = restrict_enumerative_problem(TSL,[randn(Float64,20) for i in 1:3])

visualize_with_triangles(EP;depth=3,resolution = 1000)