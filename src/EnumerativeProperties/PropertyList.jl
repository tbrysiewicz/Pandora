const system  = EnumerativeProperty{System}("system")

#We cannot define the enumerative property 'degree' since we
#   import functions from HC and Oscar called 'degree'
const enumerative_degree = EnumerativeProperty{Int}("degree")

const base_fibre = EnumerativeProperty{Fibre}("base fibre")