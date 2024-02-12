##TODO
## Implement several strategies
## Probably we should force scoring functions to be discrete x continuous
##   Otherwise, it may be too hard to have versatile code and it will be
##   too tricky for a user to start things which are not expected
## "Being Stuck", "Decay Rate", "Inflate Rate", "Taboo Tolerance", etc
##     all need to be determined by a strategy.
## There should be a TabooScore given as a score in Discrete x Cts
##    If there is no discrete part, we call the tabooscore Tolerant


export
	optimize_enumerative,
	OptimizerData,
	real_sampler,
	make_better,
	optimize_real_solns,
	optimize_reals_generic,
	default_data

## OptimizerData is the struct that contains all the information about the next sampling step.

mutable struct Strategies
	Shotgun::Bool
	Careful::Bool #Maybe it should be Tuple{Bool,Union{Float64,Nothing}} for (Bool, Radius). Think whether it's better to add radius as a part of this attribute.
	Reset::Bool
	Reveries::Bool #Maybe it should be Tuple{Bool,Union{Float64,Nothing}} for (Bool, WeightVector).
	Ambitious::Bool
end

mutable struct OptimizerData
	## This holds inside of it the record holders
	##  and information about the proportion of TABOO runs on the
	##  previous runs. This will inform the optimizer to change
	##  how it samples
	RecordFibre::Tuple{Result,Vector{Float64}}
	Record::Any
	TabooScore::Float64
	StuckScore::Int64
	PreviousFibre::Tuple{Result,Vector{Float64}}
	Radius::Float64
	WeightVector::Vector{Float64}
	Strategy::Strategies	
	#There should be a strategy flag
		#strategy: careful (push through valleys)
		#		   long-shots (additionally take large radius in a second bucket)
		#		   optimistic (increase weight of previous improvement direction)
		#Meta Strategy:
		#		   lot's of seeds
		#
end
### TODO: The only functions that might need to incorporate the Strategy attribute are real_sampler,
### 		make_better,optimizer_data_updater and optimize. Need to update them to incorporate the strategies. 



## Function for creating a sampling function that fits in with the given OptimizerData and bucket_size.

# Probably have to incorporate careful and ambitious into this. Maybe Reveries also.


#=function real_sampler(EP::EnumerativeProblem, OD::OptimizerData) ##This should also depend on OptimizerData
	k = nparameters(EP.F)		#The function being used is "nparameters" from HomotopyContinuation, not "n_parameters" defined earlier in EnumerativeProblem.jl. 
								#Hence not fitting the naming style with other functions here.
	direction = (OD.RecordFibre[2]-OD.PreviousFibre[2])
	if OD.Strategy.Careful[1] 
		println("Old Radius: ",OD.Strategy.Careful[2])
		if OD.TabooScore>0.7 #If the taboo score is too high, reduce radius
			OD.Strategy.Careful[2]=OD.Strategy.Careful[2]*0.9
		elseif OD.TabooScore<0.2 #If the taboo score is too low, increase radius
			OD.Strategy.Careful[2]=OD.Strategy.Careful[2]*1.1
		end
		println("Radius: ",OD.Strategy.Careful[2])
		function sampler1(n)
			v = [OD.RecordFibre[2]] #Always include current fibre
			for i in 1:n
				push!(v,OD.Strategy.Careful[2]*randn(Float64,k)+OD.RecordFibre[2])
			end
			if OD.Strategy.Ambitious
				for i in 1:2
					push!(v,OD.RecordFibre[2]+direction*OD.Strategy.Careful[2]*i^2) #and the fibres which continues the direction of the last move
				end
			end
			return(v)
		end
	sampler = sampler1
	else
		function sampler2(n)
			v = [OD.RecordFibre[2]]
			for i in 1:n
				push!(v,OD.Strategy.Careful[2]*randn(Float64,k)+OD.RecordFibre[2])
			end
			if OD.Strategy.Ambitious
				for i in 1:2
					push!(v,OD.RecordFibre[2]+direction*i^2)
				end
			end
			return(v)
		end
	sampler = sampler2
	end
return(sampler)
end
=#

#=
function real_sampler(EP::EnumerativeProblem, OD::OptimizerData) ##This should also depend on OptimizerData
	k = nparameters(EP.F)		#The function being used is "nparameters" from HomotopyContinuation, not "n_parameters" defined earlier in EnumerativeProblem.jl. 
								#Hence not fitting the naming style with other functions here.
	direction = (OD.RecordFibre[2]-OD.PreviousFibre[2])
	 
		println("Old Radius: ",OD.Radius)
		if OD.TabooScore>0.7 #If the taboo score is too high, reduce radius
			OD.Radius=OD.Radius*0.9
		elseif OD.TabooScore<0.2 #If the taboo score is too low, increase radius
			OD.Radius=OD.Radius*1.1
		end
		println("Radius: ",OD.Radius)
		function sampler1(n)
			v = [OD.RecordFibre[2]] #Always include current fibre
			for i in 1:n
				push!(v,OD.Radius*randn(Float64,k)+OD.RecordFibre[2])
			end
			if OD.Strategy.Ambitious
				for i in 1:2
					push!(v,OD.RecordFibre[2]+direction*OD.Radius*i^2) #and the fibres which continues the direction of the last move
				end
			end
			return(v)
		end
	return(sampler1)
end
=#

function real_sampler(EP::EnumerativeProblem, OD::OptimizerData) ##This should also depend on OptimizerData
	k = nparameters(EP.F)		#The function being used is "nparameters" from HomotopyContinuation, not "n_parameters" defined earlier in EnumerativeProblem.jl. 
								#Hence not fitting the naming style with other functions here.
	direction = (OD.RecordFibre[2]-OD.PreviousFibre[2])
	if OD.Strategy.Careful 
		println("Old Radius: ",OD.Radius)
		radius_updater(OD)
		println("Radius: ",OD.Radius)
		function sampler1(n)
			v = [OD.RecordFibre[2]] #Always include current fibre
			for i in 1:n
				push!(v,OD.Radius*randn(Float64,k)+OD.RecordFibre[2])
			end
			if OD.Strategy.Ambitious
				for i in 1:2
					push!(v,OD.RecordFibre[2]+direction*OD.Radius*i^2) #and the fibres which continues the direction of the last move for Ambitious.
				end
			end
			return(v)
		end
	sampler = sampler1
	else
		function sampler2(n)
			v = [OD.RecordFibre[2]]
			for i in 1:n
				push!(v,OD.Radius*randn(Float64,k)+OD.RecordFibre[2])
			end
			if OD.Strategy.Ambitious
				for i in 1:2
					push!(v,OD.RecordFibre[2]+direction*i^2)
				end
			end
			return(v)
		end
	sampler = sampler2
	end
return(sampler)
end

##The function for updating radius according to tabooscore, when the strategy careful is true.
function radius_updater(OD::OptimizerData)
	if OD.TabooScore>0.7 #If the taboo score is too high, reduce radius
		OD.Radius=OD.Radius*0.9
	elseif OD.TabooScore<0.2 #If the taboo score is too low, increase radius
		OD.Radius=OD.Radius*1.1
	end
end

## The function for tracking stuck scores, i.e., for tracking how long the code has been in the same 
## discriminant chamber, making very little progress.

function last_score_progress(new_rec,old_rec)
	rec_improvement = -(last(new_rec[1])-last(old_rec[1]))
	if rec_improvement<0 #This means it is larger because of another coord
		println("Progress:             ",:Inf)
		return(true)
	end
	movement = norm(new_rec[2][2]-old_rec[2][2])
	progress_score=rec_improvement/movement
	println("Progress:             ",progress_score)
	if progress_score>0.1 ###What should this number be?
		return(true)
	else
		return(false)
	end
end


## Functions for tracking tabooscores: i.e., for tracking the cases when the new samples go out of the discriminant chamber.

function first_score_taboo(sol::Tuple{Result,Vector{Float64}},SC::Score,k)
	if (SC.ScoreFunction(sol))[1]<k
		return(true)
	else
		return(false)
	end
end

function first_score_taboo_proportion(sols::Vector{Tuple{Result,Vector{Float64}}},SC::Score,OD::OptimizerData)
	taboo_counter = 0
	for sol in sols
		if first_score_taboo(sol,SC,OD.Record[1])
			taboo_counter=taboo_counter+1
		end
	end
	return(1.0*taboo_counter/length(sols))
end

## Functions for improving the sampled parameters, and updating the OptimizerData accordingly.

function make_better(EP::EnumerativeProblem,
					OD::OptimizerData,
					SC::Score;
					TS=first_score_taboo_proportion,
					PROG = last_score_progress, bucket_size=100)
	new_sampler=real_sampler(EP,OD)
	sols = solve_over_params(EP,new_sampler(bucket_size))
	#Everything below needs to be its own updating procedure
	# like, UpdateOptimizer(OD,sols)
	#This should make things easier to work with when we implement strategies
	#We could even reintroduce the sampler as an attribute of the optimizer?
	optimizer_data_updater(OD,SC,sols,TS = TS, PROG = PROG)
end

#=
Creating an update_optimizer function called optimizer_data_updater as mentioned in between the code of make_better function.
=#

#=
function optimizer_data_updater(OD::OptimizerData, SC:: Score, sols; PROG = last_score_progress, TS=first_score_taboo_proportion)
	(record,record_fibre) = max_score(sols,SC)   
	if PROG((record,record_fibre),(OD.Record,OD.RecordFibre)) && OD.Radius>0.001
		OD.StuckScore=0
	else
		OD.StuckScore=OD.StuckScore+1
	end
	OD.TabooScore = TS(sols,SC,OD)
	OD.PreviousFibre = OD.RecordFibre
	OD.Record=record
	OD.RecordFibre=record_fibre
	println("Record:", OD.Record)
	println("Taboo Score:",OD.TabooScore)
	println("Stuck Score:", OD.StuckScore)
	if OD.StuckScore>=100
		println("It seems we are stuck....going into chaos mode")
		OD.Radius=1000
		OD.StuckScore=0
	end
end
=#

function optimizer_data_updater(OD::OptimizerData, SC:: Score, sols; PROG = last_score_progress, TS=first_score_taboo_proportion)
	(record,record_fibre) = max_score(sols,SC)   
	if PROG((record,record_fibre),(OD.Record,OD.RecordFibre)) && OD.Radius>0.001
		OD.StuckScore=0
	else
		OD.StuckScore=OD.StuckScore+1
	end
	OD.TabooScore = TS(sols,SC,OD)
	OD.PreviousFibre = OD.RecordFibre
	OD.Record=record
	OD.RecordFibre=record_fibre
	println("Record:", OD.Record)
	println("Taboo Score:",OD.TabooScore)
	println("Stuck Score:", OD.StuckScore)
	if OD.StuckScore>=100
		println("It seems we are stuck....going into chaos mode")
		OD.Radius=1000
		OD.StuckScore=0
	end
end

## Finally, the general optimize function.
function optimize_enumerative(E::EnumerativeProblem, SC::Score, N;bucket_size=100)
	##First, do a really random brute force search to find a good starting point
	OD = default_data(E, SC)
	for i in 1:N
		make_better(E,OD,SC;bucket_size=bucket_size)
	end
	return(OD)
end

#=Example

AL = AreaLengthSystem()
Degree(AL)
SC = RealScoreSpace
OD = optimize(AL,SC,100)

All Real solutions!
[-1291.7598740797775, -1405.7987828197163, -1220.0662595045399, -59.1110622202469, -112.15102526093892, -371.3459096104677, -71.03932987299308, -84.27322083494848, -341.1868353357009, -1244.3188147721949]
=#



## Defining a few functions below to make it easier to check for better optimization strategies, for the case of 
##real solutions of a polynomial system. Can be changed/removed later if needed. -- Deepak

function optimize_real_solns(F::System, N; bucket_size = 100)
    SC = RealScoreSpace
    E = EnumerativeProblem(F)
    optimize_enumerative(E, SC, N, bucket_size = bucket_size)
end


function optimize_reals_generic(degree_of_poly;n_iterations=100)
	no_of_variables=degree_of_poly+1
	@var a[1:no_of_variables], x
	F = System([a[no_of_variables]+sum([a[i]*x^i for i in 1:degree_of_poly])],variables=[x],parameters=vec(a))
	optimize_real_solns(F, n_iterations)
end

## Defining an example data to make it easier to check if the code is working.
function default_data(E::EnumerativeProblem,SC::Score)
	k = nparameters(E.F)		#nparameters from HomotopyContinuation, not from Pandora. Hence the stylistic choice.
	println("Nparameters:",k)
	starting_sample_size = 1
	sols = solve_over_params(E,[randn(Float64,k) for i in 1:starting_sample_size])
	(record,record_fibre) = max_score(sols,SC)
	taboo_score = 0.5
	stuck_score = 0
	previous_fibre = record_fibre
	radius = 1
	weightvector = randn(Float64, k)
	#For defining strategy:
	shotgun = false
	careful = true
	reset = false
	reveries = false
	ambitious = true
	strategy = Strategies(shotgun,careful,reset,reveries,ambitious)
	OD = OptimizerData(record_fibre,record,taboo_score,stuck_score,previous_fibre,radius,weightvector, strategy)
	return(OD)
end




