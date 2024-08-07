export
    AltBurmester





function NewAltBurmester(M,N)
    j=M+N; #Number of points being interpolated


    @var G1r,G2r,G1i,G2i,z1r,z2r,z1i,z2i, dx[1:j],dy[1:j], p[1:j]

    #initializing the location vectors and orientation angle vectors, lowercase and uppercase pairs indicating 
    #pairs of isotropic coordinates. tr and ti are the real and imaginary parts of the orientation angles

    d=[] #location of coupler point - represented in the dx and dy variables.
    D=[] #conjugate location vector
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


    for i in 1:M-2
        push!(eq5,t[i]*z1+d[i]-g1)
        push!(eq6,T[i]*Z1+D[i]-G1)
        push!(eq8,t[i]*z2+d[i]-g2)
        push!(eq9,T[i]*Z2+D[i]-G2)
    end

	if M!=1
    for i in M:j
        push!(eq5,t[i]*z1+d[i-1]-g1)
        push!(eq6,T[i]*Z1+D[i-1]-G1)
        push!(eq8,t[i]*z2+d[i-1]-g2)
        push!(eq9,T[i]*Z2+D[i-1]-G2)
    end
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

function AltBurmester(M,N)
    j=M+N;


    @var G1r,G2r,G1i,G2i,z1r,z2r,z1i,z2i, dx[1:j],dy[1:j], theta[1:j]

    #initializing the location vectors and orientation angle vectors, lowercase and uppercase pairs indicating 
    #pairs of isotropic coordinates. tr and ti are the real and imaginary parts of the orientation angles

    d=[]
    D=[]
    t=[]
    T=[]
    tr=[]
    ti=[]
    


    for i in 1:j-1
        push!(d,(dx[i+1]-dx[1])+(dy[i+1]-dy[1])*im);
        push!(D,(dx[i+1]-dx[1])-(dy[i+1]-dy[1])*im);
    end


    #converting the Cartesian coordinates to isotropic form
    g1=G1r+im*G1i;
    g2=G2r+im*G2i;
    z1=z1r+im*z1i;
    z2=z2r+im*z2i;


    G1=G1r-im*G1i;
    G2=G2r-im*G2i;
    Z1=z1r-im*z1i;
    Z2=z2r-im*z2i;
    
    
    

    for i in 1:j
        push!(tr,1-theta[i]^2)
        push!(ti,2*theta[i])
    end
    

    for i in 1:j
        push!(t,(tr[i]+im*ti[i])/(1+theta[i]^2)) 
        push!(T,(tr[i]-im*ti[i])/(1+theta[i]^2))  
    end


    #setting equation (2), (5), (6), (7) (8), (9) (10) in 
    #the paper


    eq2 = Vector{Expression}(undef,j)
    for i in 1:j
        eq2[i] = tr[i]^2+ti[i]^2-(1+theta[i]^2)^2
    end
    
    eq5 = Vector{Expression}(undef,j)
    eq6 = Vector{Expression}(undef,j)
    eq8 = Vector{Expression}(undef,j)
    eq9 = Vector{Expression}(undef,j)


    eq5[1]=z1-g1;
    eq6[1]=Z1-G1;
    eq8[1]=z2-g2;
    eq9[1]=Z2-G2;


    for i in 1:M-1
        eq5[i+1]=t[i]*z1+d[i]-g1;
        eq6[i+1]=T[i]*Z1+D[i]-G1;
        eq8[i+1]=t[i]*z2+d[i]-g2;
        eq9[i+1]=T[i]*Z2+D[i]-G2;
    end


    for i in M:j
        eq5[i]=t[i]*z1+d[i-1]-g1;
        eq6[i]=T[i]*Z1+D[i-1]-G1;
        eq8[i]=t[i]*z2+d[i-1]-g2;
        eq9[i]=T[i]*Z2+D[i-1]-G2;
    end


    eq7=Vector{Expression}(undef,j-1)
    eq10=Vector{Expression}(undef,j-1)
    for i in 1:(j-1)
        eq7[i]=eq5[i+1]*eq6[i+1]-eq5[1]*eq6[1]
        eq10[i]=eq8[i+1]*eq9[i+1]-eq8[1]*eq9[1]
    end
    
    polys = []
    #=
    for i in 1:j
        append!(polys,eq2[i])
    end
    =#
    
    for i in 1:j-1
        append!(polys,eq7[i])
        append!(polys,eq10[i])
    end
    
    params = Vector{Variable}(undef,2*j+M)


    for i in 1:j
        params[i]=dx[i];
        params[j+i]=dy[i]
    end
    #tr[1:M] and ti[1:M] represent cos(theta) and sin(theta) respectively. they are parameters
    #in this case since theta is dependent on the poses
    for i in 1:M
        params[i+2*j]=theta[i]
        
    end




    vars=Vector{Variable}(undef,8+N)


    vars[1]=G1r;vars[2]=G2r;vars[3]=G1i;vars[4]=G2i;
    vars[5]=z1r;vars[6]=z2r;vars[7]=z1i;vars[8]=z2i;
    for i in 1:N
        vars[i+8]=theta[i+M]
    end


    F=System(polys,variables=vars,parameters=params)
    E = EnumerativeProblem(F)
    degree(E;Method="Monodromy")
    return(E)
    
end
function AltBurmester(M,N)
	j=M+N;

	mypi=3.14152697;

	@var G1r,G2r,G1i,G2i,z1r,z2r,z1i,z2i,tr[1:N],ti[1:N],theta[1:M],dx[1:j],dy[1:j]

	d=Vector{Expression}(undef,j-1)
	D=Vector{Expression}(undef,j-1)
	t=Vector{Expression}(undef,M-1)
	T=Vector{Expression}(undef,M-1)

	for i in 1:M-1
		t[i]=cos((theta[i+1]-theta[1])*mypi/180)+im*sin(( theta[i+1]-theta[1])*mypi/180);
		T[i]=cos((theta[i+1]-theta[1])*mypi/180)-im*sin(( theta[i+1]-theta[1])*mypi/180);
	end

	for i in 1:j-1
		d[i]=((dx[i+1]-dx[1])+(dy[i+1]-dy[1])*im);
		D[i]=((dx[i+1]-dx[1])-(dy[i+1]-dy[1])*im);
	end

	#converting the Cartesian coordinates to isotropic form
	g1=G1r+im*G1i;
	g2=G2r+im*G2i;
	z1=z1r+im*z1i;
	z2=z2r+im*z2i;

	G1=G1r-im*G1i;
	G2=G2r-im*G2i;
	Z1=z1r-im*z1i;
	Z2=z2r-im*z2i;

	t_var=Vector{Expression}(undef,N)
	T_var=Vector{Expression}(undef,N)
	for i in 1:N
		t_var[i]=tr[i]+im*ti[i];
		T_var[i]=tr[i]-im*ti[i];
	end

	#setting equation (2), (5), (6), (7) (8), (9) (10) in #the paper

	eq2 = Vector{Expression}(undef,N)
	for i in 1:N
		eq2[i] = t_var[i]*T_var[i]-1;
	end
	
	eq5 = Vector{Expression}(undef,j)
	eq6 = Vector{Expression}(undef,j)
	eq8 = Vector{Expression}(undef,j)
	eq9 = Vector{Expression}(undef,j)

	eq5[1]=z1-g1;
	eq6[1]=Z1-G1;
	eq8[1]=z2-g2;
	eq9[1]=Z2-G2;

	for i in 1:M-1
		eq5[i+1]=t[i]*z1+d[i]-g1;
		eq6[i+1]=T[i]*Z1+D[i]-G1;
		eq8[i+1]=t[i]*z2+d[i]-g2;
		eq9[i+1]=T[i]*Z2+D[i]-G2;
	end

	for i in 1:N
		eq5[M+i]=t_var[i]*z1+d[i+M-1]-g1;
		eq6[M+i]=T_var[i]*Z1+D[i+M-1]-G1;
		eq8[M+i]=t_var[i]*z2+d[i+M-1]-g2;
		eq9[M+i]=T_var[i]*Z2+D[i+M-1]-G2;
	end

	eq7=Vector{Expression}(undef,j-1)
	eq10=Vector{Expression}(undef,j-1)
	for i in 1:(j-1)
		eq7[i]=eq5[i+1]*eq6[i+1]-eq5[1]*eq6[1]
		eq10[i]=eq8[i+1]*eq9[i+1]-eq8[1]*eq9[1]
	end
	
	polys = []
	for i in 1:N
		append!(polys,eq2[i])
	end
	
	for i in 1:j-1
		append!(polys,eq7[i])
		append!(polys,eq10[i])
	end
	
	params = Vector{Variable}(undef,2*j+M)

	for i in 1:j
		params[i]=dx[i];
		params[j+i]=dy[i]
	end

	for i in 1:M
		params[i+2*j]=theta[i]
		
	end

	vars=Vector{Variable}(undef,8+2*N)

	vars[1]=G1r;vars[2]=G2r;vars[3]=G1i;vars[4]=G2i;
	vars[5]=z1r;vars[6]=z2r;vars[7]=z1i;vars[8]=z2i;
	for i in 1:N
		vars[i+8]=tr[i]
		vars[i+8+N]=ti[i]
	end

	F=System(polys,variables=vars,parameters=params)
	E = EnumerativeProblem(F)
    degree(E;Method="Monodromy")
    return(E)
end
#=
E=AltBurmester(M,N)
SC = RealScoreSpace

OptimizerStarters = [optimize_enumerative(E,SC,5;bucket_size=100) for i in 1:10]
OD = optimize_enumerative(E,SC,5;bucket_size=500)

howeverlongyouwant = 10
for i in 1:howeverlongyouwant
    make_better(E,OD,SC)
end


=#



#=
#The Alt Burmester problem of going through M+N points with prescribed poses at the first M of them
function AltBurmester(M,N)
	j=M+N;

	mypi=3.14152697;

	@var G1r,G2r,G1i,G2i,z1r,z2r,z1i,z2i,tr[1:N],ti[1:N],dx[1:j],dy[1:j]#,theta[1:M]

	d=Vector{Expression}(undef,j-1)
	D=Vector{Expression}(undef,j-1)
	#t=Vector{Expression}(undef,M-1)
	#T=Vector{Expression}(undef,M-1)

#=
	for i in 1:M-1
		t[i]=cos((theta[i+1]-theta[1])*mypi/180)+im*sin(( theta[i+1]-theta[1])*mypi/180);
		T[i]=cos((theta[i+1]-theta[1])*mypi/180)-im*sin(( theta[i+1]-theta[1])*mypi/180);
	end
=#

	for i in 1:j-1
		d[i]=((dx[i+1]-dx[1])+(dy[i+1]-dy[1])*im);
		D[i]=((dx[i+1]-dx[1])-(dy[i+1]-dy[1])*im);
	end

	#converting the Cartesian coordinates to isotropic form
	g1=G1r+im*G1i;
	g2=G2r+im*G2i;
	z1=z1r+im*z1i;
	z2=z2r+im*z2i;

	G1=G1r-im*G1i;
	G2=G2r-im*G2i;
	Z1=z1r-im*z1i;
	Z2=z2r-im*z2i;

    @var tr[1:M+N]
    @var ti[1:M+N]

	t_var=Vector{Expression}(undef,M+N)
	T_var=Vector{Expression}(undef,M+N)
	for i in 1:N+M
		t_var[i]=tr[i]+im*ti[i];
		T_var[i]=tr[i]-im*ti[i];
	end

	#setting equation (2), (5), (6), (7) (8), (9) (10) in #the paper

	eq2 = Vector{Expression}(undef,N)
	for i in M+1:N
		eq2[i] = t_var[i]*T_var[i]-1;
	end
	
	eq5 = Vector{Expression}(undef,j)
	eq6 = Vector{Expression}(undef,j)
	eq8 = Vector{Expression}(undef,j)
	eq9 = Vector{Expression}(undef,j)

	eq5[1]=z1-g1;
	eq6[1]=Z1-G1;
	eq8[1]=z2-g2;
	eq9[1]=Z2-G2;

	for i in 1:M+N
		eq5[i+1]=t_var[i]*z1+d[i]-g1;
		eq6[i+1]=T_var[i]*Z1+D[i]-G1;
		eq8[i+1]=t_var[i]*z2+d[i]-g2;
		eq9[i+1]=T_var[i]*Z2+D[i]-G2;
	end

	for i in 1:N
		eq5[M+i]=t_var[i]*z1+d[i+M-1]-g1;
		eq6[M+i]=T_var[i]*Z1+D[i+M-1]-G1;
		eq8[M+i]=t_var[i]*z2+d[i+M-1]-g2;
		eq9[M+i]=T_var[i]*Z2+D[i+M-1]-G2;
	end

	eq7=Vector{Expression}(undef,j-1)
	eq10=Vector{Expression}(undef,j-1)
	for i in 1:(j-1)
		eq7[i]=eq5[i+1]*eq6[i+1]-eq5[1]*eq6[1]
		eq10[i]=eq8[i+1]*eq9[i+1]-eq8[1]*eq9[1]
	end
	
	polys = []
	for i in 1:N
		append!(polys,eq2[i])
	end


	
	for i in 1:j-1
		append!(polys,eq7[i])
		append!(polys,eq10[i])
	end
	
	params = Vector{Variable}(undef,2*j+M)

    params = Vector{Variable}([])

	for i in 1:j
        push!(params,dx[i])
        push!(params,dy[i])
	end
#=
	for i in 1:M
		params[i+2*j]=theta[i]
	end
    =#
    for i in 1:M-1
        push!(params,t[i])
        push!(params,T[i])
    end

	vars=Vector{Variable}(undef,8+2*N)

	vars[1]=G1r;vars[2]=G2r;vars[3]=G1i;vars[4]=G2i;
	vars[5]=z1r;vars[6]=z2r;vars[7]=z1i;vars[8]=z2i;

	for i in 1:N
		vars[i+8]=tr[i]
		vars[i+8+N]=ti[i]
	end

	F=System(polys,variables=vars,parameters=params)
	
end
    =#