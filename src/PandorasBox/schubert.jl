export LinesMeetingFourLinesInP3

function LinesMeetingFourLinesInP3()
    @var p[1:2,1:4]
    @var q[1:4,1:2,1:4]

    Eqs = []
    for i in 1:4
        M = [
            p[1,1] p[1,2] p[1,3] p[1,4]
            p[2,1] p[2,2] p[2,3] p[2,4]
            q[i,1,1] q[i,1,2] q[i,1,3] q[i,1,4]
            q[i,2,1] q[i,2,2] q[i,2,3] q[i,2,4]
        ]
        push!(Eqs, det(M))
    end

    @var b[1:2,1:2]

    Eqs = evaluate.(
        Eqs,
        Ref(vec(p) => [
            1, 0, b[1,1], b[1,2],
            0, 1, b[2,1], b[2,2]
        ])
    )

    Params = vec(q)

    F = System(Eqs, variables=vec(b), parameters=Params)
    E = EnumerativeProblem(F)
end