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
	optimize,
	OptimizerData,
	real_sampler,
	make_better

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
	#There should be a strategy flag
		#strategy: careful (push through valleys)
		#		   long-shots (additionally take large radius in a second bucket)
		#		   optimistic (increase weight of previous improvement direction)
		#Meta Strategy:
		#		   lot's of seeds
		#
end

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
	function sampler(n)
		v = [OD.RecordFibre[2]] #Always include current fibre
		for i in 1:n
			push!(v,OD.Radius*randn(Float64,k)+OD.RecordFibre[2])
		end
		for i in 1:2
			push!(v,OD.RecordFibre[2]+direction*OD.Radius*i^2) #and the fibres which continues the direction of the last move
		end
		return(v)
	end
	return(sampler)
end

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
	if OD.StuckScore>100
		println("It seems we are stuck....going into chaos mode")
		OD.Radius=1000
		OD.StuckScore=0
	end
end


function optimize(E::EnumerativeProblem, SC::Score, N;bucket_size=100)
	##First, do a really random brute force search to find a good starting point
	k = nparameters(E.F)		#nparameters from HomotopyContinuation, not from Pandora. Hence the stylistic choice.
	println("Nparameters:",k)
	starting_sample_size = 1
	sols = solve_over_params(E,[randn(Float64,k) for i in 1:starting_sample_size])
	(record,record_fibre) = max_score(sols,SC)
	OD = OptimizerData(record_fibre,record,0.5,0,record_fibre,10000)
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

function first_score_taboo_proportion(sols::Vector{Tuple{Result,Vector{Float64}}},SC::Score,OD::OptimizerData)
	taboo_counter = 0
	for sol in sols
		if first_score_taboo(sol,SC,OD.Record[1])
			taboo_counter=taboo_counter+1
		end
	end
	return(1.0*taboo_counter/length(sols))
end

function first_score_taboo(sol::Tuple{Result,Vector{Float64}},SC::Score,k)
	if (SC.ScoreFunction(sol))[1]<k
		return(true)
	else
		return(false)
	end
end
