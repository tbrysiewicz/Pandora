export
    TwentySevenLines,
    TangentCircles,
    SymmetricTwentySevenLines




@doc raw"""
    TwentySevenLines()

 Returns the enumerative problem of finding the twenty-seven lines on a cubic surface in C^3.
 # Examples
 ```jldoctest

 julia> T = TwentySevenLines();

 julia> degree(T)
 Populating a base fibre of the enumerative problem
 Tracking 45 paths... 100%|████████████████████████████████████████| Time: 0:00:01
   # paths tracked:                  45
   # non-singular solutions (real):  27 (0)
   # singular endpoints (real):      0 (0)
   # total solutions (real):         27 (0)
 27
 
 ```
 """
function TwentySevenLines()
    @var x,y,z
    @var a[1:4,1:4,1:4]
    terms = []
    for i in 0:3
        for j in 0:3
            for k in 0:3
                if i+j+k<=3
                    push!(terms,[i,j,k])
                end
            end
        end
    end
    f = sum([a[c[1]+1,c[2]+1,c[3]+1]*x^c[1]*y^c[2]*z^c[3] for c in terms])
    Params = [a[c[1]+1,c[2]+1,c[3]+1] for c in terms]

    @var t,b[1:2],c[1:2]

    lx = t
    ly = b[1]*t+b[2]
    lz = c[1]*t+c[2]

    g = HomotopyContinuation.subs(f,[x,y,z]=>[lx,ly,lz])
    Eqs = HomotopyContinuation.coefficients(g,[t])

    F = System(Eqs,variables=[b[1],b[2],c[1],c[2]],parameters=Params)
    E = EnumerativeProblem(F)
end

function QuinticThreefolds()
    @var x[1:5]
    @var a[1:2],b[1:2],c[1:2],d[1:2],e[1:2]
    @var a50000,a41000,a32000,a31100,a22100,a21110,a11111
    f1 = [x[i]^5 for i in 1:5]
    f2 = [x[i[1]]^4*x[i[2]] for i in combinations(1:5,2)]
    f3 = [x[i[1]]^3*x[i[2]]^2 for i in combinations(1:5,2)]
    f4 = [x[i[1]]^3*x[i[2]]*x[i[3]] for i in combinations(1:5,3)]
    f5 = [x[i[1]]^2*x[i[2]]^2*x[i[3]] for i in combinations(1:5,3)]
    f6 = [x[i[1]]^2*x[i[2]]*x[i[3]]*x[i[4]] for i in combinations(1:5,4)]
    f7 = [x[1]*x[2]*x[3]*x[4]*x[5]]
    fparts = [f1,f2,f3,f4,f5,f6,f7]
    params = [a50000,a41000,a32000,a31100,a22100,a21110,a11111]
    F = sum([sum(a) for a in fparts.*params])


    @var t,a[1:2],b[1:2],c[1:2],d[1:2]

    lx = a[1]*t+a[2] #represent lines in this form t=>[lx,ly,lz,lw]
    ly = b[1]*t+b[2]
    lz = c[1]*t+c[2]
    lw = d[1]*t+d[2]
    lv = e[1]*t+e[2]
    vars = [a[1],b[1],c[1],d[1],a[2],b[2],c[2],d[2],e[1],e[2]]
    #Gr(2,5) has dimension 2*3 = 6

    g = HomotopyContinuation.subs(F,x=>[lx,ly,lz,lw,lv])
    Eqs = HomotopyContinuation.coefficients(g,[t])
    for i in 1:4
        push!(Eqs,sum(randn(Float64,10).*vars)-1.0)
    end

    Fsystem = System(Eqs,variables=vars,parameters=params)
    E = EnumerativeProblem(Fsystem)
end

function SymmetricTwentySevenLines()
    @var x,y,z,w
    @var a3000,a2100,a1110 #symmetric cubics form a 3-dimensional space
    f = a3000*(x^3+y^3+z^3+w^3) + 
         a2100*(x^2*y+x*y^2+x^2*z+x*z^2+x^2*w+w^2*x+y^2*z+z^2*y+y^2*w+w^2*y+z^2*w+w^2*z) + 
             a1110*(x*y*z+x*y*w+x*z*w+y*z*w) #they look like this
    #change_of_coords = [sum(randn(Float64,4).*[x,y,z,w]) for i in 1:4] #this is a random real change of coords so lines are not generically hiding at infinity
    #f = subs(f,[x,y,z,w]=>change_of_coords)
    f = subs(f,a3000=>1.0)
    #Params = [a3000,a2100,a1110]
    Params = [a2100,a1110]
    @var t,a[1:2],b[1:2],c[1:2],d[1:2]

    lx = a[1]*t+a[2] #represent lines in this form t=>[lx,ly,lz,lw]
    ly = b[1]*t+b[2]
    lz = c[1]*t+c[2]
    lw = d[1]*t+d[2]
    vars = [a[1],b[1],c[1],d[1],a[2],b[2],c[2],d[2]]

    g = HomotopyContinuation.subs(f,[x,y,z,w]=>[lx,ly,lz,lw])
    Eqs = HomotopyContinuation.coefficients(g,[t])
    for i in 1:4
        push!(Eqs,sum(randn(Float64,8).*vars)-1.0)
    end

    F = System(Eqs,variables=[a[1],b[1],c[1],d[1],a[2],b[2],c[2],d[2]],parameters=Params)
    E = EnumerativeProblem(F)
end

##Code taken from https://mathrepo.mis.mpg.de/circlesTangentConics/

function TangentCircles()
    @var s, t, r # Circle centered at (s, t) with radius r
    @var u[1:3], v[1:3] # Three points of tangency (u[1], v[1]), ..., (u[3], v[3])
    @var a[1:6], b[1:6], c[1:6] # Three fixed conics with coefficients defined by a[1], ..., a[6]; ... ; c[1], ..., c[6]
    
    #=
    Set up the polynomial system defining tritangent circles
    =#
    
    f_1 = (u[1] - s)^2 + (v[1] - t)^2 - r # the point (u[1], v[1]) lies on the circle
    
    f_2 = (u[2] - s)^2 + (v[2] - t)^2 - r # the point (u[2], v[2]) lies on the circle
    
    f_3 = (u[3] - s)^2 + (v[3] - t)^2 - r # the point (u[3], v[3]) lies on the circle
    
    f_4 = a[1]*u[1]^2 + a[2]*u[1]*v[1] + a[3]*v[1]^2 + a[4]*u[1] + a[5]*v[1] + a[6] # the point (u[1], v[1]) lies on the conic defined by coefficients a[1:6]
    
    f_5 = b[1]*u[2]^2 + b[2]*u[2]*v[2] + b[3]*v[2]^2 + b[4]*u[2] + b[5]*v[2] + b[6] # the point (u[2], v[2]) lies on the conic defined by coefficients b[1:6]
    
    f_6 = c[1]*u[3]^2 + c[2]*u[3]*v[3] + c[3]*v[3]^2 + c[4]*u[3] + c[5]*v[3] + c[6] # the point (u[3], v[3]) lies on the conic defined by coefficients c[1:6]
    
    f_7 = det([differentiate(f_1, [u[1], v[1]]) differentiate(f_4, [u[1], v[1]])]) # Circle and conic are tangent at (u[1], v[1])
    
    f_8 = det([differentiate(f_2, [u[2], v[2]]) differentiate(f_5, [u[2], v[2]])]) # Circle and conic are tangent at (u[2], v[2])
    
    f_9 = det([differentiate(f_3, [u[3], v[3]]) differentiate(f_6, [u[3], v[3]])]) # Circle and conic are tangent at (u[3], v[3])
    
    # Polynomial system in variables  [u[1], v[1], u[2], v[2], u[3], v[3], s, t, r] and parameters [a[1], ..., a[6], b[1], ..., b[6], c[1], ..., c[6]]
    
    F = System([f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9], variables = [u[1], v[1], u[2], v[2], u[3], v[3], s, t, r], 
    parameters = [a[1], a[2], a[3], a[4], a[5], a[6], b[1], b[2], b[3], b[4], b[5], b[6], c[1], c[2], c[3], c[4], c[5], c[6]]);
    
    return(EnumerativeProblem(F))
    end
    