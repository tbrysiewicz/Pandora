export SpaceConics


#See https://academicweb.nd.edu/~jhauenst/preprints/hsAlphaCertified.pdf
function SpaceConics(k::Int, n::Int)
    @assert k + 2n == 8 "The space of conics in P^3 has dimension 8, so need k + 2n == 8."

    @var A B C
    @var c11 c22 c12 c13 c23

    # Point parameters: p[:,i] is the i-th point in P^3.
    @var p[1:4, 1:n]

    # Line parameters: the j-th line is spanned by u[:,j] and v[:,j].
    @var u[1:4, 1:k], v[1:4, 1:k]

    variables = [A, B, C, c11, c22, c12, c13, c23]
    parameters = vcat(vec(p), vec(u), vec(v))

    # We work on the affine chart where the plane has equation
    #
    #     x4 = A*x1 + B*x2 + C*x3
    #
    # and the conic inside this plane is given by a quadratic in
    # the coordinates x1,x2,x3. We normalize the x3^2 coefficient to 1.
    q(x1, x2, x3) =
        c11*x1^2 +
        c22*x2^2 +
        x3^2 +
        c12*x1*x2 +
        c13*x1*x3 +
        c23*x2*x3

    h(x1, x2, x3, x4) = x4 - A*x1 - B*x2 - C*x3

    eqs = []

    # Point conditions:
    # A point must lie on the plane, and then satisfy the conic equation.
    for i in 1:n
        push!(eqs, h(p[1,i], p[2,i], p[3,i], p[4,i]))
        push!(eqs, q(p[1,i], p[2,i], p[3,i]))
    end

    # Line conditions:
    # The line L = span(u,v) meets the plane at the point
    #
    #     h(v)u - h(u)v.
    #
    # The conic meets L iff this intersection point satisfies q = 0.
    for j in 1:k
        hu = h(u[1,j], u[2,j], u[3,j], u[4,j])
        hv = h(v[1,j], v[2,j], v[3,j], v[4,j])

        z1 = hv*u[1,j] - hu*v[1,j]
        z2 = hv*u[2,j] - hu*v[2,j]
        z3 = hv*u[3,j] - hu*v[3,j]

        push!(eqs, q(z1, z2, z3))
    end

    return EnumerativeProblem(
        System(eqs, variables = variables, parameters = parameters)
    )
end

#=
space_conics(8, 0)  # 92 
space_conics(6, 1)
space_conics(4, 2)
space_conics(2, 3)
space_conics(0, 4)
=#