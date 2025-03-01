export 
    TwentySevenLines,
    TangentCircles
    
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

    g = subs(f,[x,y,z]=>[lx,ly,lz])
    Eqs = coefficients(g,[t])

    F = System(Eqs,variables=[b[1],b[2],c[1],c[2]],parameters=Params)
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
    