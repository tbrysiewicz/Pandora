
#Function : find_real_points
#Input: V = variety, N = number of real points to find   
#Output: Nrealpts = list containing N real points from V
function find_N_real_points(V::Variety,N::Int; limit = 100)
	D = ambient_dimension(V) 
	d = D -Pandora.dim(V) #dimension of linear subspace is codim of variety
	W = witness_set(V) 
	realpts = [] #initalizing list of real points 
	i = 0 #initalizing counter

	while i < limit 

		L = rand_subspace(D, dim=d, real = true, affine = true)
		newW = witness_set(W,L) #witness set of W using a random linear subspace
		realpts = vcat(realpts, HomotopyContinuation.real_solutions(results(newW))) #list of real points off of NewW
		length(realpts)>= N && break #break if we have found at least N real points 
		i = i+1
	end 
	Nrealpts = vcat(realpts[1:N]) #collecting N real points off of witness set
	return(Nrealpts)
end 
