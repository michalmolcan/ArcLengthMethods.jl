using Plots, BenchmarkTools, Revise
using NLsolve, LinearAlgebra, ForwardDiff

includet("../src/ArcLengthMethod.jl")
using .ArcLengthMethod

function fint(a)
    θ₀ = pi/3
    return [(1/sqrt(1-2*a[1]*sin(θ₀) + a[1]^2) -1)*(sin(θ₀) - a[1])]
end

function demo()
    fext = [1.0]
    Δl = 5e-4
    u0 = [1e-6]

    λ0=1e-2
    iterations = 10
    ftol = 1e-8

    qs = []

    ndof = length(u0)
    
    qs = []

    R = similar(u0)
    
    function evalfun!(R,u,λ)
        R[1:length(u)] = fint(u)-λ*fext
    end

    result = nlsolve((R,a) -> evalfun!(R,a,λ0),u0,method=:newton)
    u = result.zero
    push!(qs,[u;λ0])

    result = nlsolve((R,a) -> evalfun!(R,a,λ0+Δl),u,method=:newton)
    u = result.zero
    push!(qs,[u;λ0+Δl])

    λ = last(qs)[2]
    pl = scatter([last(qs)[1]],[last(qs)[2]],markershape=:auto,legend=:topleft,label="u₀")
    ylabel!("load factor λ")
    xlabel!("displacement u")

    # predictor
    qs2 = copy(qs)
    # qs2[2][1:end-1] = 5*qs2[2][1:end-1] 
    qs2[2][end] = 1.1*qs2[2][end] 

    Δq = diff(qs2)[1]
    Δq *= Δl/norm(Δq)

    Δu = Δq[1:end-1]
    Δλ = last(Δq)


    iteration = 0
    R = similar(u).+1

    Kt = ForwardDiff.jacobian(fint, u) 

    scatter!([last(qs)[1]]+Δu,[last(qs)[2]+Δλ],markershape=:auto, label="iteration 0")

    converged = false
    while iteration < iterations-1
        evalfun!(R,u + Δu,λ + Δλ)
        converged = (norm(R)/ndof < ftol)
        if converged; break; end

        iteration += 1
        Δu,Δλ = mcrcorrection(Δu,Δλ,Kt,R,fext,Δl)
        scatter!([last(qs)[1]]+Δu,[last(qs)[2]+Δλ],markershape=:auto,label="iteration " * string(iteration))

        
        if isnan(Δλ); converged = false; break; end
    end

    return pl
end

demo()