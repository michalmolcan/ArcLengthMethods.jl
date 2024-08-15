using Plots, Revise
using NLsolve, LinearAlgebra, ForwardDiff

includet("../src/ArcLengthMethods.jl")
using .ArcLengthMethods

function fint(a)
    θ₀ = pi/3
    return [(1/sqrt(1-2*a[1]*sin(θ₀) + a[1]^2) -1)*(sin(θ₀) - a[1])]
end

function plotiterations(correctorstep!, methodname; ftol=1e-3, is_cylindrical = false)

    markersize = 6

    fext = [1.0]
    Δl = 2.5e-1
    u0 = [1e-6]

    iterations = 10

    qs = []

    ndof = length(u0)
    
    qs = []

    R = similar(u0)
    
    function evalfun!(R,u,λ)
        R[1:length(u)] = fint(u)-λ*fext
    end

    
    qs = arclengthmethod(fint,fext,1e-2,u0;verbose=false,method=:riks,adaptivestep=false)
    pl = plot([u[1] for u in qs],[u[2] for u in qs],ls=:auto,legend = :outertopright,label="Solution",linealpha=0.5)
    ylabel!("load factor λ")
    xlabel!("displacement u")

    stepselect = 45

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

    Δu = @view Δq[1:end-1]
    Δλ = @view Δq[end]

    push!(us,u+Δu)
    push!(λs,λ.+Δλ)
    scatter!([u+Δu],[λ.+Δλ],markershape=:auto, label="predictor",markeralpha=.4,markersize=markersize)

    iteration = 0
    bffr = (
        Δur = similar(u0),
        Δuf = similar(u0),
        Δλbar = similar(Δλ),
        α = zeros(1)
    )
    Kt = ForwardDiff.jacobian(fint, u) 

    converged = false
    while iteration < iterations-1
        evalfun!(R,u + Δu,λ .+ Δλ)
        converged = (norm(R)/ndof < ftol)
        if converged; break; end

        iteration += 1
        correctorstep!(Δu,Δλ,Kt,R,fext,Δl,bffr;cylindrical=is_cylindrical)

        push!(us,u+Δu)
        push!(λs,λ.+Δλ)
        scatter!([u+Δu],[λ.+Δλ],markershape=:auto,label="iteration " * string(iteration),markeralpha=.4,markersize=markersize)

        
        if isnan(first(Δλ)); converged = false; break; end
    end

    xlims!(minimum(us)[1]-2.5e-2,maximum(us)[1]+5e-2)
    ylims!(minimum(λs)-2.5e-2,maximum(λs)+5e-2)

    title!(methodname * " iterations, ftol = " * string(ftol))
    return pl
end

opts = (
    ftol = 1e-3,
    cylindrical = false
)

pl1 = plotiterations(rikscorrection!,"Riks method";opts...)
pl2 = plotiterations(crisfieldcorrection!, "Crisfields method";opts...)
pl3 = plotiterations(rammcorrection!, "Ramms method";opts...)
pl4 = plotiterations(mcrcorrection!, "MCR method";opts...)

savefig(pl1,"docs/files/riksmethoditerations.png")
savefig(pl2,"docs/files/crisfieldsmethoditerations.png")
savefig(pl3,"docs/files/rammsmethoditerations.png")
savefig(pl4,"docs/files/mcrmethoditerations.png")