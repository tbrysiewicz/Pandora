#=
Nodes on symmetroids in P^3.

This code constructs the enumerative problem of the nodes of a degree d
symmetroid surface in P^3.  The surface is defined as the determinant of a
linear symmetric matrix pencil

    A(x) = x_1 I_d + x_2 A_1 + x_3 A_2 + x_4 A_3,

where A_1, A_2, A_3 are general symmetric d x d matrices.  The symmetroid is

    det(A(x)) = 0.

A node of a general symmetroid occurs at a projective point x where A(x) has
corank 2.  Equivalently, x is a singular point of the determinant hypersurface,
so it is cut out by the four partial derivatives of det(A(x)).  Since these
equations are homogeneous in x, we add a random affine chart to select one
representative of each projective solution.

For real parameter values, we classify the real nodes by signature.  A real
node is called spectrahedral if the nonzero eigenvalues of A(x) are all of the
same sign, equivalently if the rank-(d-2) part of A(x) is definite up to the
projective ambiguity x ~ -x.  Thus the signature of a real solution set is

    (# real nodes, # spectrahedral nodes).

For d = 4, a general quartic symmetroid has 10 nodes.  The stratification of
real quartic spectrahedra by the location of these nodes is studied in:

    J. C. Ottem, K. Ranestad, B. Sturmfels, C. Vinzant,
    "Quartic Spectrahedra", Mathematical Programming 151 (2015), 585--612.

For d = 5, a general quintic symmetroid has 20 nodes.  The corresponding
classification by real/spectrahedral node configurations is studied in:

    T. Brysiewicz, K. Kozhasov, M. Kummer,
    "Nodes on quintic spectrahedra", Le Matematiche 76 (2021), 415--430.


=#

function nodes_on_symmetroid(d::Int)

    @var x[1:4], a[1:3,1:d,1:d]

    Id = Matrix{Int}(LinearAlgebra.I, d, d)

    A = vcat(
        [Id],
        [
            vcat([
                hcat([a[k,min(i,j),max(i,j)] for i in 1:d]...)
                for j in 1:d
            ]...)
            for k in 1:3
        ]
    )

    Ax = sum([x[i]*A[i] for i in 1:4])
    D = det(Ax)

    affine_chart = sum([randn(Float64)*x[i] for i in 1:4]) + randn(Float64)

    eqs = vcat(
        [differentiate(D, x[i]) for i in 1:4],
        [affine_chart]
    )

    symm_parameters = vcat([
        a[i,j,k]
        for i in 1:3
        for j in 1:d
        for k in j:d
    ]...)

    F = System(eqs; variables=x, parameters=symm_parameters)
    return EnumerativeProblem(F)
end

function is_spectrahedral_node(d::Int, s, p; eig_tol=1e-7)
    if !Pandora.is_real(s)
        return false
    end

    M = real(symmetroid_matrix_from_solution(d, real.(collect(s)), real.(collect(p))))
    λ = eigvals(M)

    num_positive = count(x -> x > eig_tol, λ)
    num_negative = count(x -> x < -eig_tol, λ)

    return num_positive == d - 2 || num_negative == d - 2
end

function symmetroid_signature(d::Int, sols, p; eig_tol=1e-7)
    num_real = 0
    num_spectrahedral = 0

    for s in sols
        if Pandora.is_real(s)
            num_real += 1

            if is_spectrahedral_node(d, s, p; eig_tol=eig_tol)
                num_spectrahedral += 1
            end
        end
    end

    return (num_real, num_spectrahedral)
end

function symmetroid_matrix_from_solution(d::Int, s, p)
    # s = solution coordinates for x[1:4]
    # p = parameter values ordered as:
    #     a[i,j,k] for i in 1:3 for j in 1:d for k in j:d

    xvals = collect(s)
    pvals = collect(p)

    Ms = Matrix{eltype(pvals)}[]

    push!(Ms, Matrix{eltype(pvals)}(LinearAlgebra.I, d, d))

    idx = 1
    for ell in 1:3
        M = zeros(eltype(pvals), d, d)

        for j in 1:d
            for k in j:d
                M[j,k] = pvals[idx]
                M[k,j] = pvals[idx]
                idx += 1
            end
        end

        push!(Ms, M)
    end

    return sum(xvals[i] * Ms[i] for i in 1:4)
end