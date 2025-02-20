#We cannot define the enumerative property 'degree' since we
#   import functions from HC and Oscar called 'degree'
const enumerative_degree = EnumerativeProperty{Int}("degree")

const system  = EnumerativeProperty{System}("system")

const base_fibre = EnumerativeProperty{Fibre}("base fibre")
const bkk_bound = EnumerativeProperty{Int}("bkk torus bound")
const affine_bkk_bound = EnumerativeProperty{Int}("bkk affine bound")
const degree_sequence = EnumerativeProperty{Vector{Int}}("degree sequence")
const bezout_bound = EnumerativeProperty{Int}("bezout bound")
const newton_polytopes = EnumerativeProperty{Vector{Polyhedron}}("newton polytopes")