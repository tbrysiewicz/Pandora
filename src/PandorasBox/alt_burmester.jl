#As formulated in "The complete solution of Alt-Burmester synthesis problems for four-bar linkages"
#by Brake et. al. 
#doi: https://doi.org/10.1115/1.4033251


function AltBurmester(M,N)
    #Cartesian coordinates
    @var g1[1:2] g2[1:2] z1[1:2] z2[1:2] t[1:M+N] T[1:M+N] Tbar[1:M+N] dx[1:M+N] dy[1:M+N]
    eqs = []
    #map to isotropic coordinates
    G1 = g1[1] + im*g1[2]
    G1bar = g1[1] - im*g1[2]

    G2 = g2[1] + im*g2[2]
    G2bar = g2[1] - im*g2[2]

    Z1 = z1[1] + im*z1[2]
    Z1bar = z1[1] - im*z1[2]

    Z2 = z2[1] + im*z2[2]
    Z2bar = z2[1] - im*z2[2]

    D = [(dx[j]-dx[1]) + im*(dy[j]-dy[1]) for j in 1:M+N]
    Dbar = [(dx[j]-dx[1]) - im*(dy[j]-dy[1]) for j in 1:M+N]

    # T = [(t[j] - im) / (t[j] + im) for j in 1:M+N]
    # Tbar = [(t[j] + im) / (t[j] - im) for j in 1:M+N]

    # adding equations for angles
    for j in 1:M+N
        push!(eqs, T[j]*Tbar[j] - 1)
        push!(eqs, T[j]*(t[j]+im)- (t[j]-im))
    end

    # equations for proximal links
    # since we have offset by the first precision point, the equation relating the position of the proximal links at D1 is different from the rest
    L1 = [(D[j] + T[j]*Z1 - G1) for j in 1:M+N]
    L1bar = [(Dbar[j] + Tbar[j]*Z1bar - G1bar) for j in 1:M+N]
    L2 = [(D[j] + T[j]*Z2 - G2) for j in 1:M+N]
    L2bar = [(Dbar[j] + Tbar[j]*Z2bar - G2bar) for j in 1:M+N]

    for j in 2:M+N
        push!(eqs, L1[j]*L1bar[j] - L1[1]*L1bar[1])
        push!(eqs, L2[j]*L2bar[j] - L2[1]*L2bar[1])
    end
    if M == 4 && N ==2
        F = System(eqs, parameters = vcat(dx, dy, t[1:2], t[5:6]), variables = vcat(g1, g2, z1, z2, t[3:4], T, Tbar))
    else
        F = System(eqs, parameters = vcat(dx, dy, t[1:M]), variables = vcat(g1, g2, z1, z2, t[M+1:M+N], T, Tbar))
    end

    return F
end

