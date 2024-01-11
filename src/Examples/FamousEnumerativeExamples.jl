export
    TwentySevenLines




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