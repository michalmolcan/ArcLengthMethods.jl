using Plots, BenchmarkTools, Revise
using NLsolve, LinearAlgebra, ForwardDiff

includet("../src/ArcLengthMethod.jl")
using .ArcLengthMethod

function fint(a)
    θ₀ = pi/3
    return [(1/sqrt(1-2*a[1]*sin(θ₀) + a[1]^2) -1)*(sin(θ₀) - a[1])]
end

function plotiterations(correctionstep,methodname)
    markersize = 6

    fext = [1.0]
    Δl = 2.5e-1
    u0 = [1e-6]

    λ0=.2
    iterations = 10
    ftol = 1e-3

    qs = []

    ndof = length(u0)
    
    qs = []

    R = similar(u0)
    
    function evalfun!(R,u,λ)
        R[1:length(u)] = fint(u)-λ*fext
    end

    qs = arclengthmethod(fint,fext,rikscorrection,1e-2,u0;adaptivestep=false)
    pl = plot([u[1] for u in qs],[u[2] for u in qs],ls=:auto,legend = :outertopright,label="Solution",linealpha=0.5)
    ylabel!("load factor λ")
    xlabel!("displacement u")

    stepselect = 50

    u = qs[stepselect][1:end-1]
    λ = qs[stepselect][end]

    function circleShape(h,k,r)
        θ = LinRange(0,2pi,100)
        h .+ r*sin.(θ), k .+ r*cos.(θ)
    end

    plot!(circleShape(u,λ,Δl),label="Δl",linealpha=0.5)

    scatter!([u],[λ],markershape=:auto,label="u₀",markersize=markersize)
    

    us = [u]
    λs = [λ]

    # predictor
    Δq = qs[stepselect] - qs[stepselect-1]
    Δq *= Δl/norm(Δq)

    Δu = Δq[1:end-1]
    Δλ = last(Δq)

    push!(us,u+Δu)
    push!(λs,λ+Δλ)
    scatter!([u+Δu],[λ+Δλ],markershape=:auto, label="predictor",markeralpha=.4,markersize=markersize)

    iteration = 0
    Kt = ForwardDiff.jacobian(fint, u) 

    converged = false
    while iteration < iterations-1
        evalfun!(R,u + Δu,λ + Δλ)
        converged = (norm(R)/ndof < ftol)
        if converged; break; end

        iteration += 1
        Δu,Δλ = correctionstep(Δu,Δλ,Kt,R,fext,Δl)

        push!(us,u+Δu)
        push!(λs,λ+Δλ)
        scatter!([u+Δu],[λ+Δλ],markershape=:auto,label="iteration " * string(iteration),markeralpha=.4,markersize=markersize)

        
        if isnan(Δλ); converged = false; break; end
    end

    xlims!(minimum(us)[1]-2.5e-2,maximum(us)[1]+5e-2)
    ylims!(minimum(λs)-2.5e-2,maximum(λs)+5e-2)

    title!(methodname * " iterations, ftol = " * string(ftol))
    return pl
end

pl1 = plotiterations(rikscorrection,"Riks method")
pl2 = plotiterations(crisfieldcorrection, "Crisfields method")
pl3 = plotiterations(rammcorrection, "Ramms method")
pl4 = plotiterations(mcrcorrection, "MCR method")

savefig(pl1,"docs/files/riksmethoditerations.png")
savefig(pl2,"docs/files/crisfieldsmethoditerations.png")
savefig(pl3,"docs/files/rammsmethoditerations.png")
savefig(pl4,"docs/files/mcrmethoditerations.png")