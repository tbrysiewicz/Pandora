
#From https://arxiv.org/pdf/2603.25003
function twisted_cubic_bisecants()
    @var m[1:4,1:4]
    @var t[1:2] s[1:2]

    # v(u) = (1,u,u^2,u^3)^T
    v(u) = [one(u), u, u^2, u^3]

    # M*v(u)
    function Mv(u, M)
        [sum(M[i,j] * v(u)[j] for j in 1:4) for i in 1:4]
    end

    # 3x3 determinant
    det3(A) = A[1,1]*(A[2,2]*A[3,3] - A[2,3]*A[3,2]) -
              A[1,2]*(A[2,1]*A[3,3] - A[2,3]*A[3,1]) +
              A[1,3]*(A[2,1]*A[3,2] - A[2,2]*A[3,1])

    vt1 = v(t[1])
    vt2 = v(t[2])
    Ms1 = Mv(s[1], m)
    Ms2 = Mv(s[2], m)

    # rows 1,2,3
    A113_s1 = [
        vt1[1] vt2[1] Ms1[1]
        vt1[2] vt2[2] Ms1[2]
        vt1[3] vt2[3] Ms1[3]
    ]
    A113_s2 = [
        vt1[1] vt2[1] Ms2[1]
        vt1[2] vt2[2] Ms2[2]
        vt1[3] vt2[3] Ms2[3]
    ]

    # rows 1,2,4
    A124_s1 = [
        vt1[1] vt2[1] Ms1[1]
        vt1[2] vt2[2] Ms1[2]
        vt1[4] vt2[4] Ms1[4]
    ]
    A124_s2 = [
        vt1[1] vt2[1] Ms2[1]
        vt1[2] vt2[2] Ms2[2]
        vt1[4] vt2[4] Ms2[4]
    ]

    f1 = det3(A113_s1)
    f2 = det3(A113_s2)
    f3 = det3(A124_s1)
    f4 = det3(A124_s2)

    F = System(
        [f1, f2, f3, f4],
        variables = [t[1], t[2], s[1], s[2]],
        parameters = vec(m),
    )

    return EnumerativeProblem(F)
end
