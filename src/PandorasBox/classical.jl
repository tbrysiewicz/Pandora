export LinesOnTwoQuadricsInP4, PlaneConicsOnQuinticThreefold


function LinesOnTwoQuadricsInP4()
    # Lines on a complete intersection of two quadrics in P^4.
    # Expected number: 16.

    @var x, y, z, w
    @var a[1:2,1:3,1:3,1:3,1:3]

    terms = []

    for i in 0:2
        for j in 0:2
            for k in 0:2
                for l in 0:2
                    if i + j + k + l <= 2
                        push!(terms, [i,j,k,l])
                    end
                end
            end
        end
    end

    f = []

    for s in 1:2
        push!(
            f,
            sum(
                a[s,c[1]+1,c[2]+1,c[3]+1,c[4]+1] *
                x^c[1] * y^c[2] * z^c[3] * w^c[4]
                for c in terms
            )
        )
    end

    Params = Variable[]

    for s in 1:2
        for c in terms
            push!(Params, a[s,c[1]+1,c[2]+1,c[3]+1,c[4]+1])
        end
    end

    @var t
    @var b[1:2], c[1:2], d[1:2]

    lx = t
    ly = b[1]*t + b[2]
    lz = c[1]*t + c[2]
    lw = d[1]*t + d[2]

    Eqs = []

    for s in 1:2
        g = evaluate(f[s], [x,y,z,w] => [lx,ly,lz,lw])
        append!(Eqs, coefficients(g, [t]))
    end

    Vars = vcat(b, c, d)

    F = System(Eqs, variables=Vars, parameters=Params)
    E = EnumerativeProblem(F)
end




function PlaneConicsOnQuinticThreefold()
    # Plane conics on a quintic threefold in P^4.
    # Expected number: 609250.
    #
    # Model:
    #   1. Choose an affine chart on Gr(3,5), i.e. a plane P^2 in P^4.
    #   2. Choose a conic q(u,v)=0 in that plane.
    #   3. Restrict the quintic F to the plane.
    #   4. Require F|_{P^2} to be divisible by q:
    #
    #          F|_{P^2} = q * h
    #
    #      where h is a cubic in u,v.
    #
    # Variables:
    #   plane variables: 6
    #   conic variables: 5, after fixing coefficient of v^2 to 1
    #   cubic multiplier variables: 10
    # Total variables: 21
    #
    # Equations:
    #   coefficient comparison in degree <= 5 in u,v gives 21 equations.

    @var x, y, z, w
    @var a[1:6,1:6,1:6,1:6]

    # General affine quintic in A^4.
    terms5 = []

    for i in 0:5
        for j in 0:5
            for k in 0:5
                for l in 0:5
                    if i + j + k + l <= 5
                        push!(terms5, [i,j,k,l])
                    end
                end
            end
        end
    end

    F5 = sum(
        a[c[1]+1,c[2]+1,c[3]+1,c[4]+1] *
        x^c[1] * y^c[2] * z^c[3] * w^c[4]
        for c in terms5
    )

    Params = Variable[]

    for c in terms5
        push!(Params, a[c[1]+1,c[2]+1,c[3]+1,c[4]+1])
    end

    # Affine coordinates on the moving plane.
    @var u, v

    # Plane in A^4:
    #   x = u
    #   y = v
    #   z = b11*u + b12*v + b13
    #   w = b21*u + b22*v + b23
    @var b[1:2,1:3]

    lx = u
    ly = v
    lz = b[1,1]*u + b[1,2]*v + b[1,3]
    lw = b[2,1]*u + b[2,2]*v + b[2,3]

    restricted_F = evaluate(F5, [x,y,z,w] => [lx,ly,lz,lw])

    # General affine conic in the plane, in the chart where coeff(v^2)=1.
    @var q[1:5]

    conic =
        q[1] +
        q[2]*u +
        q[3]*v +
        q[4]*u^2 +
        q[5]*u*v +
        v^2

    # General affine cubic multiplier h(u,v).
    cubic_terms = []

    for i in 0:3
        for j in 0:3
            if i + j <= 3
                push!(cubic_terms, [i,j])
            end
        end
    end

    @var h[1:10]

    cubic = sum(
        h[r] * u^cubic_terms[r][1] * v^cubic_terms[r][2]
        for r in 1:length(cubic_terms)
    )

    diff = restricted_F - conic*cubic

    Eqs = coefficients(diff, [u,v])

    Vars = vcat(vec(b), q, h)

    F = System(Eqs, variables=Vars, parameters=Params)
    E = EnumerativeProblem(F; monodromy=true)
end