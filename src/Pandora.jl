module Pandora

using Pkg 					#To remove eventually

using HomotopyContinuation	
using LinearAlgebra
using Combinatorics
using Oscar					


import HomotopyContinuation.solve

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
          ..................
				TCB
				");

		print("Version")
		printstyled(" $VERSION_NUMBER ", color = :green)
end




include("CoreObjects/EnumerativeProblem.jl")
include("CoreMethods/MovingInParameterSpace.jl")

include("Enumeration/DegreeBounds.jl")

include("GaloisGroups/GaloisSampling.jl")
include("GaloisGroups/GaloisGroups.jl")

include("CoreObjects/Variety.jl")

include("AlgebraicMatroid/AlgebraicMatroid.jl")

include("Exploration/Explore.jl")
include("Exploration/FibreFunctions.jl")
include("Exploration/Samplers.jl")

include("Examples/BasicEnumerativeExamples.jl")
include("Examples/FamousEnumerativeExamples.jl")
include("Examples/BasicVarietyExamples.jl")

include("Statistics/Tally.jl")


include("Optimization/CoreObjects.jl")
include("Optimization/Optimizer.jl")



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
