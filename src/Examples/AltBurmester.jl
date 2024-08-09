export
    AltBurmester





function AltBurmester(M,N)
    j=M+N; #Number of points being interpolated


    @var G1r,G2r,G1i,G2i,z1r,z2r,z1i,z2i, dx[1:j],dy[1:j], p[1:j]

    #initializing the location vectors and orientation angle vectors, lowercase and uppercase pairs indicating 
    #pairs of isotropic coordinates. tr and ti are the real and imaginary parts of the orientation angles

    d=[dx[1]-dx[1]+(dy[1]-dy[1])*im] #location of coupler point - represented in the dx and dy variables.
    D=[dx[1]-dx[1]-(dy[1]-dy[1])*im] #conjugate location vector
    t=[] #orientation angles
    T=[] #conjugate orientation angles
    tr=[] #real part of orientation angle
    ti=[] #imaginary part of orientation angle
    

	#Construct the coupler point locations as functions of dx and dy
    for i in 1:j-1
        push!(d,(dx[i+1]-dx[1])+(dy[i+1]-dy[1])*im);
        push!(D,(dx[i+1]-dx[1])-(dy[i+1]-dy[1])*im);
    end


    #g1,g2 are the pinned point vectors at first interpolation
    g1=G1r+im*G1i;
    g2=G2r+im*G2i;
    z1=z1r+im*z1i; #z1/z2 are circle point vectors at first interpolation (others are functions of these and the theta angles)
    z2=z2r+im*z2i;


    G1=G1r-im*G1i;
    G2=G2r-im*G2i;
    Z1=z1r-im*z1i;
    Z2=z2r-im*z2i;
    
    
    

    for i in 1:j #tr+iti = re^itheta but tr and ti are parametrized by a rational parametrization of the circle (denominators cleared)
        push!(tr,1-p[i]^2)
        push!(ti,2*p[i])
    end
    

    for i in 1:j
        push!(t,(tr[i]+im*ti[i])/(1+p[i]^2)) 
        push!(T,(tr[i]-im*ti[i])/(1+p[i]^2))  
    end

    #setting equation (2), (5), (6), (7) (8), (9) (10) in 
    #the paper


    eq2 = []
    for i in 1:j
		push!(eq2,tr[i]^2+ti[i]^2-(1+p[i]^2)^2)
    end
    

    eq5=[z1-g1];
    eq6=[Z1-G1];
    eq8=[z2-g2];
    eq9=[Z2-G2];


    for i in 1:j
        push!(eq5,t[i]*z1+d[i]-g1)
        push!(eq6,T[i]*Z1+D[i]-G1)
        push!(eq8,t[i]*z2+d[i]-g2)
        push!(eq9,T[i]*Z2+D[i]-G2)
    end

	for e in eq5
		println(e)
	end


    eq7=[]
    eq10=[]
    for i in 1:(j-1)
        push!(eq7,eq5[i+1]*eq6[i+1]-eq5[1]*eq6[1])
        push!(eq10,eq8[i+1]*eq9[i+1]-eq8[1]*eq9[1])
    end
    
    polys = []
    
    for i in 1:j-1
        push!(polys,eq7[i])
        push!(polys,eq10[i])
    end
    
    params = Vector{Variable}([])


    for i in 1:j
		push!(params,dx[i])
	end
	for i in 1:j
		push!(params,dy[i])
    end

    #tr[1:M] and ti[1:M] represent cos(theta) and sin(theta) respectively. they are parameters
    #in this case since theta is dependent on the poses
    for i in 1:M
        push!(params,p[i])
    end




    vars=Vector{Variable}([G1r,G2r,G1i,G2i,z1r,z2r,z1i,z2i])


    for i in 1:N
        push!(vars,p[i+M])
    end


    F=System(polys,variables=vars,parameters=params)
    E = EnumerativeProblem(F)
    degree(E;Method="Monodromy")
    return(E)
end


