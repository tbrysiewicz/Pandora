
export KuramotoModel


################################################
########### Problem Formulation ################
################################################


function KuramotoModel(n)
    @var w[1:(n-1)], s[1:n], c[1:n]
    equations = []
    for i in 1:(n-1)
        sum = 0
        for j in 1:n
            sum += (s[i]*c[j] - s[j]*c[i])
        end
        f_1 = w[i] - (1/n)*sum
        f_2 = c[i]^2 + s[i]^2 - 1
        f_1 = subs(f_1, [s[n], c[n]]=>[0, 1])
        f_2 = subs(f_2, [s[n], c[n]]=>[0, 1])
        f_1 == 0 || push!(equations, f_1)
        f_2 == 0 || push!(equations, f_2)
    end
    F = System(equations, variables = [s[1:n-1]..., c[1:n-1]...], parameters = [w[1:n-1]...])

    return EnumerativeProblem(F)
end



################################################
####### Specialized Fibre Functions ############
################################################







################################################
####### (s,p) Visualization Functions ##########
################################################






