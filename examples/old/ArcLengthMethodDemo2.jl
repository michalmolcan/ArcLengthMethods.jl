using NLsolve, LinearAlgebra, Plots, ForwardDiff, BenchmarkTools

function newtonraphsonapproach(nlfun!,λs,u)
    us = []

    for λ in λs
        result = nlsolve((R,a) -> nlfun!(R,a,λ),u,method =:newton)
        u = result.zero
        push!(us,u)
    end

    λs,us
end

function arclengthcontinuation(nlfun!,fext,Δl)

    λ0 = 0.01
    u0 = [0,1e-6]

    function paranlfun!(R,q,(q₀,Δl))
        u = q[1:end-1]
        λ = last(q)
    
        nlfun!(R,u,λ)
        R[3] = transpose(q-q₀)*(q-q₀) - Δl^2
    end

    qs = []

    result = nlsolve((R,a) -> nlfun!(R,a,λ0),u0,method =:newton)
    u = result.zero
    push!(qs,[u;λ0])


    result = nlsolve((R,a) -> nlfun!(R,a,λ0+Δl),u,method =:newton)
    u = result.zero
    push!(qs,[u;λ0+Δl])

    idx = 1
    while last(qs)[end] < 1
        q = secantPredictor(last(qs,2),Δl)
        result = nlsolve((R,q) -> paranlfun!(R,q,(last(qs),Δl)),q,method =:newton)
        push!(qs,result.zero)
        idx += 1
        if idx > 1000; break; end
    end

    return qs
end


function arclengthmethod(fint,fext,correctorstep,Δl;ftol=1e-8,iterations=100,steps=1000)

    λ0 = 0.01
    u0 = [0,1e-6]
    ndof = length(u0)

    qs = []

    R = similar(u0).+1

    function nlfun!(R,u,λ)
        R[1:2] = fint(u)-λ*fext
    end

    result = nlsolve((R,a) -> nlfun!(R,a,λ0),u0,method=:newton)
    u = result.zero
    push!(qs,[u;λ0])

    result = nlsolve((R,a) -> nlfun!(R,a,λ0+Δl),u,method=:newton)
    u = result.zero
    push!(qs,[u;λ0+Δl])
    
    
    step = 1
    while last(qs)[end] < 1
        
        u = last(qs)[1:end-1]
        λ = last(qs)[end]

        q = secantpredictor(last(qs,2),Δl)
        Δu = q[1:end-1]-u
        Δλ = last(q)-λ

        Kt = ForwardDiff.jacobian(fint, u)
        
        iteration = 0
        R = similar(u).+1

        converged = false
        while iteration < iterations-1
            nlfun!(R,u + Δu,λ + Δλ)
            converged = (norm(R)/ndof < ftol)
            if converged; break; end

            iteration += 1
            Δu,Δλ = correctorstep(Δu,Δλ,Kt,R,fext,Δl)
            
            if isnan(Δλ); converged = false; break; end
        end

        if converged 
            u = u + Δu
            λ = λ + Δλ
            
            push!(qs,[u;λ])
            println("Step ",step, " converged at ",iteration, " iterations, λ=",λ)

            Δl = updatearclength(Δl,iteration)

            step +=1
        else
            Δl = Δl/2
            println("Step diverged, Δl reduced to ",Δl)
        end
        if step > steps; break; end
    end

    qs
end

function secantpredictor(qs,Δl)
    Δq = diff(qs)[1]
    return last(qs) + Δl/norm(Δq)*Δq
end


function updatearclength(Δl,iterations)
    Δl *= 4/max(iterations,1)
end

"""
# Riks corrector step
[1] Ferreira2005
"""
function rikscorrection(Δu,Δλ,Kt,R,fext,Δl)
    φ = 1

    A = [Kt                 -fext;
        transpose(2*Δu)     2*Δλ*φ^2*transpose(fext)*fext]
    δq = -A\[R;transpose(Δu)*Δu + Δλ^2*φ^2*transpose(fext)*fext - Δl^2]

    Δu = Δu + δq[1:end-1]
    Δλ = Δλ + last(δq) 

    Δu,Δλ
end


"""
# Crisfields corrector step
[1] Ferreira2005
"""
function crisfieldcorrection(Δu,Δλ,Kt,R,fext,Δl)
    δū = -Kt\R
    δut = Kt\fext

    φ = 1

    a = transpose(δut)*δut + φ^2*transpose(fext)*fext
    b = 2*transpose(δut)*(Δu + δū) + 2*Δλ*φ^2*transpose(fext)*fext
    c = transpose(Δu + δū)*(Δu + δū) - Δl^2 + Δλ^2*φ^2*transpose(fext)*fext

    D = b^2 - 4*a*c
    @assert D >= 0 "Discriminant must not be negative!"

    δλ = (-b .+ [-1;1]*sqrt(D))/2/a
    Δu1 = Δu + δū + δλ[1]*δut
    Δu2 = Δu + δū + δλ[2]*δut
    
    if transpose(Δu1)*Δu/Δl^2 >= transpose(Δu2)*Δu/Δl^2
        Δλ = Δλ + δλ[1]
        Δu = Δu1
    else
        Δλ = Δλ + δλ[2]
        Δu = Δu2
    end

    Δu,Δλ
end

"""
# Ramm corrector step
[1] Fafard1993
"""
function rammcorrection(Δu,Δλ,Kt,R,fext,Δl)
    Δur = -Kt\R
    Δuf = Kt\fext

    Δλ2 = -(transpose(Δu)*Δur)/(transpose(Δu)*Δuf)
    
    Δu += Δur + Δλ2*Δuf
    Δλ += Δλ2 

    Δu,Δλ
end



function fint(a)
    E = 210e9
    A = pi*(5e-3)^2
    l₀ = 4
    L₀ = 1
    θ₀ = pi/3

    β = E*A/l₀
    k = E*A/L₀

    w = β/k

    B(a₁,θ₀) = 1-2*a₁*sin(θ₀) + a₁^2

    F₁(a₁,a₂) = (1/sqrt(B(a₁,θ₀))-1)*(sin(θ₀)-a₁) - w*(a₂ - a₁)
    F₂(a₁,a₂) = w*(a₂ - a₁)

    return [F₁(a[1],a[2]);F₂(a[1],a[2])]
end

function equilibriumstate!(R,u,λ,fext)
    R[1:2] = fint(u)-λ*fext
end


function demo1()

    # problem definition
    fext = [0, 1] 
    λs = 0.01:0.01:0.6

    nlfun(R,a,λ) = equilibriumstate!(R,a,λ,fext)
    u0 = [0,1e-6]

    # Newthon-Raphson method
    λs,us = newtonraphsonapproach(nlfun,λs,u0)
    pl = plot([u[2] for u in us],λs,label="NR approach",legend=:bottomright)
    ylabel!("load factor λ")
    xlabel!("displacement u")

    println("NR done")

    Δl = 5.0e-2

    # Arc length continuation
    qs = arclengthcontinuation(nlfun,fext,Δl)
    plot!([u[2] for u in qs],[u[3] for u in qs],ls=:auto,label="Parameterization")

    

    # Riks arc length method
    println("Riks arc length method")
    qs = arclengthmethod(fint,fext,rikscorrection,Δl)
    println("Riks method done")
    plot!([u[2] for u in qs],[u[3] for u in qs],ls=:auto,label="Riks")

    # Crisfields arc length method
    println("Crisfields arc length method")
    qs = arclengthmethod(fint,fext,crisfieldcorrection,Δl)
    println("Crisfields method done")
    plot!([u[2] for u in qs],[u[3] for u in qs],ls=:auto,label="Crisfield")

    # Ramm arc length method
    println("Ramms arc length method")
    qs = arclengthmethod(fint,fext,rammcorrection,Δl)
    println("Ramms method done")
    plot!([u[2] for u in qs],[u[3] for u in qs],ls=:auto,label="Ramm")

    # Modified Crisfield-Ramm arc length method
    # qs = arclengthmethod(mcrcorrector)
    # plot!([u[2] for u in qs],[u[3] for u in qs],ls=:auto,label="MCR")

    return pl
end

demo1()

# # Benchmark
# @btime arclengthcontinuation($nlfun,$fext,$Δl)
# @btime arclengthmethod($nlfun,$fext,$rikscorrectorstep,$Δl)
# @btime arclengthmethod($nlfun,$fext,$crisfieldcorrectorstep,$Δl)
# @btime arclengthmethod($nlfun,$fext,$rammcorrectorstep,$Δl)