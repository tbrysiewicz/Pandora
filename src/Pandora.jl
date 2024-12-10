module Pandora

using Pkg 					#remove eventually


#LinearAlgebra functions used
using LinearAlgebra: I, norm

#HC functions used (in the code)
using HomotopyContinuation:System,  @var,  expressions,  monodromy_solve, Result, TrackerOptions
using HomotopyContinuation:subs, coefficients

#HC functions which we extend to different types
import HomotopyContinuation:variables, parameters, system, solutions, solve, is_real

#HC functions we want the user to have access to
export System, @var, solve, is_real

#Oscar (GAP) functions used (in code)
using Oscar: perm, PermGroupElem, symmetric_group, sub, gens, PermGroup, cperm
using Oscar: is_transitive, is_primitive, describe, orbits, minimal_block_reps
using Oscar: cycles, orbit, on_sets, order
#Oscar (GAP) functions which we extend
import Oscar: degree


#Oscar (GAP) functions we want users to have access to
export gens, sub, symmetric_group, perm, cperm, cycles
export is_transitive, is_primitive, describe, orbits, minimal_block_reps
export orbit,  transitivity, order

#Base
import Base: convert




#Exports of Pandora types
export
    AbstractEnumerativeProblem,
    EnumerativeProblem,
    Fibre

#Exports of Pandora functions
export
    system,
    n_parameters,
    n_polynomials,
    base_solutions,
    solutions,
    base_parameters,
    parameters,
    base_fibre,
    ambient_dimension,
    degree,
    solve!,
    monodromy_solve!,
    is_populated,
    populate!,
    monodromy_homomorphism,
    monodromy_sample,
    monodromy_group,
    galois_group,
    is_decomposable,
    is_irreducible,
    components,
    data,
    tracker_options,
    restrict,
    orbit_block_partition,
    orbit_block_relabel!,
    transitivity_basis,
    monodromy_dictionary,
    monodromy_basis,
    is_score

export
    TwentySevenLines

export 
    Sampler,
    EllipseSampler,
    sample,
    initialize_real_sampler,
    dietmaier,
    initialize_real_optimizer


#using HomotopyContinuation	
#using LinearAlgebra
#using Combinatorics
#using Oscar			
#using Plots
#using ProgressBars
#using Clustering
#using ImplicitPlots



#import HomotopyContinuation.solve
#import HomotopyContinuation.degree
#import HomotopyContinuation.dim
#import Oscar.dim
#import Oscar.degree

#export order

const PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
const VERSION_NUMBER = VersionNumber(PROJECT_TOML["version"])


function __init__()
	print(raw"
        0ooo000oo oo oo
     0oo00o0o0o  0o0 o8oo
     0000o0o000o00o0000000
     ooo 0oo 0 oo00ooo00000
          oo00  /o o  ooooo
            \\\//  /////
               \\\////
                 |||/\
                 |||\/
                 |||||
           .....//||||\....
         _______/PANDO(RA)\_____
          ...................
				TCB
				");

		print("Version")
		printstyled(" $VERSION_NUMBER ", color = :green)
end

include("Constants/TypeAliases.jl")
include("Constants/Warnings.jl")

include("DependencyConversions/HC_Conversions.jl")
include("DependencyConversions/Julia_Conversions.jl")
include("DependencyConversions/Oscar_GAP_Conversions.jl")

include("EnumerativeProblems/AbstractEnumerativeProblem.jl")


include("EnumerativeProblems/EnumerativeProblem.jl")
include("EnumerativeProblems/EPGetters.jl")
include("EnumerativeProblems/EPSetters.jl")
include("EnumerativeProblems/EPSolving.jl")
include("EnumerativeProblems/EPWrappers.jl")

include("Fibres/Fibre.jl")
include("Fibres/FibreFunctions.jl")

include("Samplers/Sampler.jl")


include("Monodromy/MonodromyGroups.jl")
include("Monodromy/MonodromyInterpretations.jl")
include("Monodromy/MonodromyLabeling.jl")
include("Monodromy/MonodromyBreakup.jl")


include("Examples/NamedExamples.jl")


include("SolutionFunctions/SolutionFunctions.jl")

include("Optimization/Optimizer.jl")
include("Optimization/OptimizationGetters.jl")
include("Optimization/DietmaierOptimization.jl")

#include("Varieties/Variety.jl")
#include("Fibres/Scores.jl")

#This needs to be fixed in several ways.
#include("AlgebraicMatroid/AlgebraicMatroid.jl")

#include("CoreMethods/MovingInParameterSpace.jl")


#include("EnumerativeProblems/EPGetters.jl")

#include("Enumeration/DegreeBounds.jl")

#include("Examples/BasicEnumerativeExamples.jl")
#include("Examples/FamousEnumerativeExamples.jl")
#include("Examples/BasicVarietyExamples.jl")
#include("Examples/MatroidRealizations.jl")
#include("Examples/CrossRatioMap.jl")
#include("Visualization/MatroidVisualization.jl")
#include("Examples/AltBurmester.jl")
#include("Examples/CCEquations.jl")


#include("Exploration/Explore.jl")
#include("Exploration/FibreFunctions.jl")
#include("Exploration/Samplers.jl")

#include("Fibres/FibreChecks.jl")
#include("Fibres/FibreFunctions.jl")
#include("Fibres/RealityScores.jl")

#include("GaloisGroups/GaloisSampling.jl")
#include("GaloisGroups/GaloisGroups.jl")
#include("GaloisGroups/CoordinateSymmetryGroup.jl")

#include("Optimization/CoreObjects.jl")
#include("Optimization/Optimizer.jl")
#include("Optimization/OptimizerBeta.jl")
#include("Optimization/OptimizerUpdaters.jl")


#include("Statistics/Tally.jl")




#include("Visualization/Visualize.jl")
#include("Visualization/Triangulate.jl")
#include("Visualization/Visualize_Discriminant.jl")

#include("Sampling/Sampling.jl")

#include("Automation/Automate.jl")

end # module Pandora


#=
Doc template

@doc raw"""
    b(n)

 Returns
 # Examples
 ```jldoctest

 ```
 """
=#
