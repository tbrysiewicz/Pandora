module Pandora

using Pkg 					#To remove eventually

using HomotopyContinuation	
using LinearAlgebra
using Combinatorics
using Oscar			
using Plots
using ProgressBars
using Clustering
using ImplicitPlots



import HomotopyContinuation.solve
import HomotopyContinuation.degree
import HomotopyContinuation.dim
import Oscar.dim
import Oscar.degree

export order

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

include("CoreObjects/Varieties/Variety.jl")
include("CoreObjects/EnumerativeProblems/EnumerativeProblem.jl")
#include("Fibres/Scores.jl")

#This needs to be fixed in several ways.
#include("AlgebraicMatroid/AlgebraicMatroid.jl")

include("CoreMethods/MovingInParameterSpace.jl")


include("CoreObjects/EnumerativeProblems/EPGetters.jl")

include("Enumeration/DegreeBounds.jl")

include("Examples/BasicEnumerativeExamples.jl")
include("Examples/FamousEnumerativeExamples.jl")
include("Examples/BasicVarietyExamples.jl")
#include("Examples/MatroidRealizations.jl")
#include("Visualization/MatroidVisualization.jl")
include("Examples/AltBurmester.jl")


include("Exploration/Explore.jl")
include("Exploration/FibreFunctions.jl")
include("Exploration/Samplers.jl")

include("Fibres/FibreChecks.jl")
include("Fibres/FibreFunctions.jl")
include("Fibres/RealityScores.jl")

include("GaloisGroups/GaloisSampling.jl")
include("GaloisGroups/GaloisGroups.jl")
include("GaloisGroups/CoordinateSymmetryGroup.jl")

#include("Optimization/CoreObjects.jl")
#include("Optimization/Optimizer.jl")
include("Optimization/OptimizerBeta.jl")


include("Statistics/Tally.jl")




include("Visualization/Visualize.jl")
include("Visualization/Triangulate.jl")




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
