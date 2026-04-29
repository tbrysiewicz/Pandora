export BarthSextic, SphereProjection


function BarthSextic()
    @var w, x, y, z
    f = 4*(w^2*x^2 - y^2)*(w^2*y^2 - z^2)*(w^2*z^2 - x^2) - (1 + 2*w)*(x^2 + y^2 + z^2 - 1^2)^2
    golden_ratio = (1 + sqrt(5))/2
    f = subs(f, w=>golden_ratio)
    F = System([f], variables = [z], parameters = [x,y])

    return EnumerativeProblem(F)
end

function SphereProjection()
    @var x,y,z
    @var a,b
    F = System([x^2+y^2+z^2-1,x-a,y-b],variables = [x,y,z], parameters = [a,b])
    return(EnumerativeProblem(F))
end
    
