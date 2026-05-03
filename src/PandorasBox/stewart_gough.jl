
export StewartGough

################################################
########### Problem Formulation ################
################################################


# The formulation is from the Dietmaier's 1998 paper titled :
# "The Stewart–Gough platform of general geometry can have 40 real postures" 
# In: *Advances in Robot Kinematics: Analysis and Control*. Springer, Dordrecht. https://doi.org/10.1007/978-94-015-9064-8_1.

"""
    StewartGough()

    Formulates the forward kinematics problem for the general Stewart-Gough platform.
    
**Reference:**
Dietmaier, P. (1998). The Stewart-Gough platform of general geometry can have 40 real postures. In: *Advances in Robot Kinematics: Analysis and Control*. Springer, Dordrecht. https://doi.org/10.1007/978-94-015-9064-8_1
"""
function StewartGough()
    # Define Variables (n, e1, e2)
    @var n[1:3] e1[1:3] e2[1:3]
    vars = [n; e1; e2]

    # Define Parameters
    # Define L as a vector of length 6. L[1] is unused (since L1=1).
    @var L[1:6]
    @var a21 b21 a31 a32 b31 b32
    # a and b matrices for joints 4, 5, and 6
    @var a[1:3, 1:3] b[1:3, 1:3]
    
    # Construct the parameter vector (29 total parameters)
    # We use L[2:6] to exclude the fixed length of the first leg.
    params = [L[2:6]; a21; b21; a31; a32; b31; b32; vec(a); vec(b)]

    # Define derived vector e3 = e1 x e2
    e3 = [
        e1[2]*e2[3] - e1[3]*e2[2],
        e1[3]*e2[1] - e1[1]*e2[3],
        e1[1]*e2[2] - e1[2]*e2[1]
    ]

    # Construct Equations
    eqs = []

    # Leg 2 distance restriction (Eq. 1)
    push!(eqs, sum((n + b21*e1 - [a21, 0, 0]).^2) - L[2]^2)

    # Leg 3 distance restriction (Eq. 2)
    push!(eqs, sum((n + b31*e1 + b32*e2 - [a31, a32, 0]).^2) - L[3]^2)

    # Legs 4, 5, 6 distance restrictions (Eqs. 3, 4, 5)
    for i in 1:3
        # i=1 -> Joint 4, Leg 4
        # i=2 -> Joint 5, Leg 5
        # i=3 -> Joint 6, Leg 6
        platform_joint = b[i,1]*e1 + b[i,2]*e2 + b[i,3]*e3
        base_joint = [a[i,1], a[i,2], a[i,3]]
        
        # Accessing L[4], L[5], and L[6]
        push!(eqs, sum((n + platform_joint - base_joint).^2) - L[i+3]^2)
    end

    # Geometric constraints (Eqs. 6 to 9)
    push!(eqs, sum(n.^2) - 1)      # L1 = 1 normalization
    push!(eqs, sum(e1.^2) - 1)     # e1 is unit vector
    push!(eqs, sum(e2.^2) - 1)     # e2 is unit vector
    push!(eqs, sum(e1 .* e2))      # e1 perpendicular to e2

    return EnumerativeProblem(System(eqs, variables=vars, parameters=params))
end