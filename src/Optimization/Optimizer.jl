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
	Optimize,
	OptimizerData,
	RealSampler,
	MakeBetter

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

function RealSampler(EP::EnumerativeProblem, OD::OptimizerData) ##This should also depend on optimizerdata
	k = nparameters(EP.F)
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

function LastScoreProgress(newRec,oldRec)
	RecImprovement = -(last(newRec[1])-last(oldRec[1]))
	if RecImprovement<0 #This means it is larger because of another coord
		println("Progress:             ",:Inf)
		return(true)
	end
	Movement = norm(newRec[2][2]-oldRec[2][2])
	ProgressScore=RecImprovement/Movement
	println("Progress:             ",ProgressScore)
	if ProgressScore>0.1 ###What should this number be?
		return(true)
	else
		return(false)
	end
end

function MakeBetter(EP::EnumerativeProblem,
					OD::OptimizerData,
					SC::Score;
					TS=FirstScoreTabooProportion,
					PROG = LastScoreProgress, bucketsize=100)
	newsampler=RealSampler(EP,OD)
	Sols = solve_over_params(EP,newsampler(bucketsize))
	#Everything below needs to be its own updating procedure
	# like, UpdateOptimizer(OD,Sols)
	#This should make things easier to work with when we implement strategies
	#We could even reintroduce the sampler as an attribute of the optimizer?
	(Record,RecordFibre) = MaxScore(Sols,RealScoreSpace)   
	if PROG((Record,RecordFibre),(OD.Record,OD.RecordFibre)) && OD.Radius>0.001
		OD.StuckScore=0
	else
		OD.StuckScore=OD.StuckScore+1
	end
	OD.TabooScore = TS(Sols,SC,OD)
	OD.PreviousFibre = OD.RecordFibre
	OD.Record=Record
	OD.RecordFibre=RecordFibre
	println("Record:", OD.Record)
	println("Taboo Score:",OD.TabooScore)
	println("Stuck Score:", OD.StuckScore)
	if OD.StuckScore>100
		println("It seems we are stuck....going into chaos mode")
		OD.Radius=1000
		OD.StuckScore=0
	end
end


function Optimize(E::EnumerativeProblem, SC::Score, N;bucketsize=100)
	##First, do a really random brute force search to find a good starting point
	k = nparameters(E.F)
	println("Nparameters:",k)
	StartingSampleSize = 1
	Sols = solve_over_params(E,[randn(Float64,k) for i in 1:StartingSampleSize])
	(Record,RecordFibre) = MaxScore(Sols,SC)
	OD = OptimizerData(RecordFibre,Record,0.5,0,RecordFibre,10000)
	for i in 1:N
		MakeBetter(E,OD,SC;bucketsize=bucketsize)
	end
	return(OD)
end

#=Example

AL = AreaLengthSystem()
Degree(AL)
SC = RealScoreSpace
OD = Optimize(AL,SC,100)

All Real solutions!
[-1291.7598740797775, -1405.7987828197163, -1220.0662595045399, -59.1110622202469, -112.15102526093892, -371.3459096104677, -71.03932987299308, -84.27322083494848, -341.1868353357009, -1244.3188147721949]
=#

function FirstScoreTabooProportion(Sols::Vector{Tuple{Result,Vector{Float64}}},SC::Score,OD::OptimizerData)
	taboocounter = 0
	for Sol in Sols
		if FirstScoreTaboo(Sol,SC,OD.Record[1])
			taboocounter=taboocounter+1
		end
	end
	return(1.0*taboocounter/length(Sols))
end

function FirstScoreTaboo(Sol::Tuple{Result,Vector{Float64}},SC::Score,k)
	if (SC.ScoreFunction(Sol))[1]<k
		return(true)
	else
		return(false)
	end
end
