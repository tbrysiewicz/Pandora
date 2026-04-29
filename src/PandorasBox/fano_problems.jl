
function LinesOnQuintic3Fold()
    @var x[1:5], a[1:6, 1:6, 1:6, 1:6, 1:6]

    terms = [(i, j, k, l, m) for i in 0:5, j in 0:5, k in 0:5, l in 0:5, m in 0:5 if i + j + k + l + m == 5]

    f = sum(
        a[i + 1, j + 1, k + 1, l + 1, m + 1] *
        x[1]^i * x[2]^j * x[3]^k * x[4]^l * x[5]^m
        for (i, j, k, l, m) in terms
    )

    affine_f = subs(f, x[5] => 1)

    @var t, b[1:2, 1:5]
    line = [
        b[1, 1] + t * b[2, 1],
        b[1, 2] + t * b[2, 2],
        b[1, 3] + t * b[2, 3],
        t,
    ]

    g = subs(affine_f, [x[1], x[2], x[3], x[4]] => line)

    params = [a[i + 1, j + 1, k + 1, l + 1, m + 1] for (i, j, k, l, m) in terms]
    F = System(
        coefficients(g, [t]),
        variables = [b[1, 1], b[1, 2], b[1, 3], b[2, 1], b[2, 2], b[2, 3]],
        parameters = params,
    )

    return EnumerativeProblem(F)
end

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

function SymmetricTwentySevenLines()
    @var x,y,z
    @var a[1:4,1:4,1:4,1:4]
    terms = []
    for i in 0:3
        for j in 0:3
            for k in 0:3
                if i+j+k<=3
                    push!(terms,[i,j,k, 3-i-j-k])
                end
            end
        end
    end

    f = sum([a[sort([c[1]+1,c[2]+1,c[3]+1,c[4]+1])...]*x^c[1]*y^c[2]*z^c[3] for c in terms])
    Params = unique([a[sort([c[1]+1,c[2]+1,c[3]+1,c[4]+1])...] for c in terms])

    @var t,d[1:2],b[1:2],c[1:2]

    lx = d[1]*t+d[2]
    ly = b[1]*t+b[2]
    lz = c[1]*t+c[2]

    g = subs(f,[x,y,z]=>[lx,ly,lz])
    Eqs = coefficients(g,[t])
    push!(Eqs, sum(randn(Float64,6).*vcat(d,b,c))-1)
    push!(Eqs, sum(randn(Float64,6).*vcat(d,b,c))-1)

    F = System(Eqs,variables=[b[1],b[2],c[1],c[2],d[1],d[2]],parameters=Params)
    E = EnumerativeProblem(F)
end



