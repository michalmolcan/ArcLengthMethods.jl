# function stiffnessmatrix(u)
#     E = 2.0e11
#     r = 15.0e-3
#     mu = 0.3

#     x₀ = [0;0;1000;-1000]*1e-3

#     A = pi*r^2
#     x = x₀+u;

#     L₀ = sqrt((x₀[3]-x₀[1])^2 + (x₀[4]-x₀[2])^2)
#     L = sqrt((x[3]-x[1])^2 + (x[4]-x[2])^2);

#     cs = (x[3] - x[1])/L;
#     ss = (x[4] - x[2])/L;
#     T = [cs ss  0   0; 
#         -ss cs  0   0; 
#         0   0   cs  ss; 
#         0   0   -ss cs]

#     FF = E*A/L₀*((u[3]-u[1]) + 1/2/L₀*(u[4]-u[2])^2);

#     k_E = E*A/L*[1 0 -1 0; 
#                 0 0 0 0;
#                 -1 0 1 0; 
#                 0 0 0 0];

#     k_G = FF/L*[0 0 0 0; 
#                 0 1 0 -1; 
#                 0 0 0 0; 
#                 0 -1 0 1];
    
#     return T'*(k_E + k_G)*T;
# end

# function main()

#     u₁ = [0;0;0;1]*1e-9
#     F = [0;0;0;1e8]


#     K(u) = stiffnessmatrix([0;0;0;u])[4,4]
#     f(u,λ) = K(u)*u - λ*F[4]

#     lambdas = 0.01:0.01:1
#     us = []

#     u = [u₁[4]]
#     for λ in lambdas
#         function g!(R,u)
#             R[1] = f(u[1],λ)
#         end

#         result = nlsolve(g!,u)

#         u = result.zero
#         push!(us,u[1])
#     end

#     lambdas*F[4],us
# end


# lambdas,us = main()

# plot(lambdas,us)
# xlabel!("load factor λ")
# ylabel!("displacement u")