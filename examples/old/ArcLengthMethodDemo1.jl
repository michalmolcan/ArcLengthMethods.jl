# Demonstration of arc length method, 
# code without refactoring and optimization

using NLsolve, LinearAlgebra, Plots, ForwardDiff

function newtonraphsonapproach(nlfun!,λs,u)
    us = []

    for λ in λs
        result = nlsolve((R,a) -> nlfun!(R,a,λ),u,method =:newton)
        u = result.zero
        push!(us,u)
    end

    λs,us
end

function arclengthcontinuation()
    f = [0, 1]

    λ₀ = 0.01
    u₀ = [0,1e-6]
    arc = 1e-2

    qs = []

    result = nlsolve((R,a) -> equilibriumstate!(R,a,λ₀,f),u₀,method =:newton)
    u = result.zero
    push!(qs,[u;λ₀])


    result = nlsolve((R,a) -> equilibriumstate!(R,a,λ₀+arc,f),u,method =:newton)
    u = result.zero
    push!(qs,[u;λ₀+arc])

    idx = 1
    while last(qs)[end] < 1
        q = secantPredictor(last(qs,2),arc)
        result = nlsolve((R,q) -> parameterizednlfun!(R,q,(last(qs),arc,f)),q,method =:newton)
        push!(qs,result.zero)
        idx += 1
        if idx > 1000; break; end
    end
    
    # function correctorFun(u,α,Kₜ,F,arc)
    #     rK = Kₜ(u)*u - α*F
    #     rARC = (q-qₚ)'*(q-qₚ) .- arc^2
    
    #     return [rK;rARC]
    # end

    return qs
end


function arclengthmethod(nlfun!,fext,correctorfun)

    λ₀ = 0.01
    u₀ = [0,1e-6]
    Δl = 1e-2

    qs = []

    result = nlsolve((R,a) -> nlfun!(R,a,λ₀),u₀,method =:newton)
    u = result.zero
    push!(qs,[u;λ₀])


    result = nlsolve((R,a) -> nlfun!(R,a,λ₀+Δl),u,method =:newton)
    u = result.zero
    push!(qs,[u;λ₀+Δl])
    
    R = similar(first(qs))

    step = 1
    while last(qs)[end] < 1
        p₀ = last(qs)[1:end-1]
        λ₀ = last(qs)[end]

        q = secantPredictor(last(qs,2),Δl)
        Δp = q[1:end-1]-p₀
        Δλ = last(q)-λ₀

        Kₜ = ForwardDiff.jacobian(fint, p₀)
        
        p, λ = correctorfun(nlfun!, Kₜ, p₀,Δp,λ₀,Δλ,fext,Δl)

        push!(qs,[p;λ])

        step +=1
        if step > 1000; break; end
    end

    qs
end

"""
# Corrector stet by Riks 

Ferreira2005
"""
function rikscorrector(Ψ!,Kₜ,p₀,Δp,λ₀,Δλ,fext,Δl;ftol=1e-8,iterations=100)
    R = similar(p₀).+1

    iteration = 1
    while norm(R) > ftol
        Ψ!(R,p₀ + Δp,λ₀ + Δλ)

        φ = 1

        A = [Kₜ -fext;
            transpose(2*Δp) 2*Δλ*φ^2*transpose(fext)*fext]
        δq = -A\[R;transpose(Δp)*Δp + Δλ^2*φ^2*transpose(fext)*fext - Δl^2]

        Δp = Δp + δq[1:end-1]
        Δλ = Δλ + last(δq) 
        
        iteration += 1
        if iteration > iterations; break; end
    end

    p = p₀ + Δp
    λ = λ₀ + Δλ

    return p,λ
end


"""
# Corrector step by Crisfield (Ferreira2005)

Ψ! is residual function 
"""
function crisfieldscorrector(Ψ!,Kₜ,p₀,Δp,λ₀,Δλ,fext,Δl;ftol=1e-8,iterations=100)

    maxiters = 100
    iters = 1
    R = similar(p₀).+1

    while norm(R) > 1e-3
        Ψ!(R,p₀ + Δp,λ₀ + Δλ)

        δp̄ = -Kₜ\R
        δpₜ = Kₜ\fext

        φ = 1

        a = transpose(δpₜ)*δpₜ + φ^2*transpose(fext)*fext
        b = 2*transpose(δpₜ)*(Δp + δp̄) + 2*Δλ*φ^2*transpose(fext)*fext
        c = transpose(Δp + δp̄)*(Δp + δp̄) - Δl^2 + Δλ^2*φ^2*transpose(fext)*fext

        D = b^2 - 4*a*c
        @assert D >= 0 "Discriminant must not be negative!"

        δλ = (-b .+ [-1;1]*sqrt(D))/2/a
        Δp₁ = Δp + δp̄ + δλ[1]*δpₜ
        Δp₂ = Δp + δp̄ + δλ[2]*δpₜ
        
        if transpose(Δp₁)*Δp/Δl^2 >= transpose(Δp₂)*Δp/Δl^2
            Δλ = Δλ + δλ[1]
            Δp = Δp₁
        else
            Δλ = Δλ + δλ[2]
            Δp = Δp₂
        end
        
        iters += 1
        if iters > maxiters; break; end
    end

    p = p₀ + Δp
    λ = λ₀ + Δλ

    return p,λ
end


"""
# Corrector step by Ramm (fafard1993)

Modified Riks-Wempner method (Riks1979,Wempner1971)

"""
function rammcorrector(Ψ!,Kₜ,u₀,Δu,λ₀,Δλ,fext,Δl;ftol=1e-8,iterations=100)
    iters = 1

    u = u₀ + Δu
    λ = λ₀ + Δλ

    R = similar(u).+1

    while norm(R) > 1e-3
        Ψ!(R,u,λ)

        Δur = -Kₜ\R
        Δuf = Kₜ\fext

        Δλ = -(transpose(Δu)*Δur)/(transpose(Δu)*Δuf)
        λ = λ + Δλ

        Δu = Δur + Δλ*Δuf

        u = u + Δu

        iters += 1
        if iters > iterations; break; end
    end

    return u,λ
end

"""
Modified Crisfield-Ramm method (fafard1993)

"""
function mcrcorrector(Ψ!,Kₜ,u₀,Δu,λ₀,Δλ,fext,Δl;ftol=1e-8,iterations=100)
    u = u₀ + Δu
    λ = λ₀ + Δλ

    R = similar(u).+1
    iteration = 1
    while norm(R) > ftol
        Ψ!(R,u,λ)

        Δur = -Kₜ\R
        Δuf = Kₜ\fext

        Δλ = -(transpose(Δu)*Δur)/(transpose(Δu)*Δuf)
        λ = λ + Δλ

        Δu = Δur + Δλ*Δuf

        α = Δl^2/(transpose(Δu)*Δu)
        
        u = u + α*Δu

        iteration += 1
        if iteration > iterations; break; end
    end

    return u,λ
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

function parameterizednlfun!(R,q,(q₀,Δl,fext))
    u = q[1:end-1]
    λ = last(q)

    equilibriumstate!(R,u,λ,fext)
    R[3] = transpose(q-q₀)*(q-q₀) - Δl^2
end

function secantPredictor(qs,Δl)
    Δq = diff(qs)[1]
    return last(qs) + Δl/norm(Δq)*Δq
end


function demonstration1()

    # problem definition
    fext = [0, 1] 
    λs = 0.01:0.01:0.6

    nlfun!(R,u,λ) = equilibriumstate!(R,u,λ,fext)
    u0 = [0,1e-6]

    # Newthon-Raphson method
    λs,us = newtonraphsonapproach(nlfun!,λs,u0)
    plot([u[2] for u in us],λs,label="NR approach",legend=:bottomright)
    ylabel!("load factor λ")
    xlabel!("displacement u")


    # Arc length continuation
    qs = arclengthcontinuation()
    plot!([u[2] for u in qs],[u[3] for u in qs],ls=:auto,label="Parameterization")

    # Riks arc length method
    qs = arclengthmethod(nlfun!,fext,rikscorrector)
    plot!([u[2] for u in qs],[u[3] for u in qs],ls=:auto,label="Riks")


    # Crisfields arc length method
    qs = arclengthmethod(nlfun!,fext,crisfieldscorrector)
    plot!([u[2] for u in qs],[u[3] for u in qs],ls=:auto,label="Crisfield")

    # Ramm arc length method
    qs = arclengthmethod(nlfun!,fext,rammcorrector)
    plot!([u[2] for u in qs],[u[3] for u in qs],ls=:auto,label="Ramm")

    # Modified Crisfield-Ramm arc length method
    # qs = arclengthmethod(mcrcorrector)
    # plot!([u[2] for u in qs],[u[3] for u in qs],ls=:auto,label="MCR")

end

demonstration1()




# Demonstration