module ArcLengthMethods

using NLsolve, LinearAlgebra, ForwardDiff

export arclengthmethod
export rikscorrection!,crisfieldcorrection!,rammcorrection!,mcrcorrection!

function arclengthmethod(fint,fext,Δl,u0;
    λ0=1e-2,ftol=1e-8,method=:crisfield,cylindrical=true,
    max_iters=100,steps=1000,
    verbose=false,adaptivestep=true,save_iters = false)
    
    if method == :riks
        correctorstep! = rikscorrection!
    elseif method == :crisfield
        correctorstep! = crisfieldcorrection!
    elseif method == :ramm
        correctorstep! = rammcorrection!
    elseif method == :mcr
        correctorstep! = mcrcorrection!
    else
        error("Unknown method")
    end

    if save_iters
        data = []
    end

    ndof = length(u0)
    
    qs = Vector{Float64}[]
    R = similar(u0)
    λ = zeros(1)
    Δq = zeros(ndof+1)
    Δu = @view Δq[1:end-1]
    Δλ = @view Δq[end]
    bffr = (
        Δur = similar(u0),
        Δuf = similar(u0),
        Δλbar = similar(Δλ),
        α = zeros(1)
    )

    function evalfun!(R,u,λ)
        R[1:length(u)] = fint(u) - λ.*fext
    end

    result = nlsolve((R,a) -> evalfun!(R,a,λ0),u0,method=:newton)
    u = result.zero
    push!(qs,[u;λ0])
    if save_iters
        push!(data,(1,copy(u0),[copy(λ0)],u-u0,0,0,[0]))
    end

    result = nlsolve((R,a) -> evalfun!(R,a,λ0+Δl),u,method=:newton)
    u2 = result.zero
    push!(qs,[u2;λ0+Δl])
    if save_iters
        push!(data,(1,copy(u),[copy(λ0+Δl)],u2-u,0,0,[0]))
    end

    step = 1
    while last(qs)[end] < 1
        
        u .= @view last(qs)[1:end-1]
        λ .= @view last(qs)[end]

        secantpredictor!(Δq,qs,Δl)

        if method == :riks
            Kt = ForwardDiff.jacobian(q -> fint(q[1:end]) .- last(q).*fext, [u + Δu;λ+ Δλ]) 
        else
            Kt = lu(ForwardDiff.jacobian(fint, u + Δu))
        end

        iteration = 0
        R = similar(u).+1

        converged = false
        while iteration < max_iters-1 
            evalfun!(R,u + Δu,λ + Δλ)
            if save_iters
                push!(data,(iteration+1,copy(u),copy(λ),copy(Δu),copy(Δλ),copy(Δl),copy(R)))
            end

            converged = (norm(R)/ndof < ftol)
            if converged; break; end

            iteration += 1
            correctorstep!(Δu,Δλ,Kt,R,fext,Δl,bffr;cylindrical)
            
            if isnan(Δλ[1]); converged = false; break; end
        end

        if converged 
            u = u + Δu
            λ = λ + Δλ
            
            push!(qs,[u;λ])
            if verbose; println("Step ",step, " converged at ",iteration, " iterations, λ=",λ);end

            if adaptivestep
                Δl = updatearclength(Δl,iteration)
            end

            step +=1
        else
            Δl = Δl/2
            if verbose; println("Step diverged, Δl reduced to ",Δl); end
        end
        if step > steps; break; end
    end

    if save_iters
        return qs,data
    else
        return qs
    end
end

function secantpredictor!(Δq,qs,Δl)
    @views Δq .= qs[end] - qs[end-1]
    Δq .= Δl/norm(Δq)*Δq
end


function updatearclength(Δl,iterations)
    Δl *= (4/max(iterations,1))
end

"""
# Riks corrector step
[1] Ferreira2005
"""
function rikscorrection!(Δu,Δλ,Kt,R,fext,Δl,bffr;cylindrical=false)
    if cylindrical
        φ = 0
    else
        φ = 1
    end
    
    A = [Kt                ;
        transpose(2*Δu)     2*Δλ.*φ^2*transpose(fext)*fext]
    δq = -A\[R;transpose(Δu)*Δu + Δλ.^2.0.*φ^2*transpose(fext)*fext - Δl^2]

    @views Δu .= Δu + δq[1:end-1]
    @views Δλ .= Δλ .+ δq[end]

    nothing
end


"""
# Crisfields corrector step
[1] Ferreira2005
"""
function crisfieldcorrection!(Δu,Δλ,Kt,R,fext,Δl,(;Δur,Δuf);cylindrical=false)
    Δur .= -(Kt\R)
    Δuf .= Kt\fext

    η = 1

    if cylindrical
        φ = 0
    else
        φ = 1
    end

    a = transpose(Δuf)*Δuf .+ φ^2*transpose(fext)*fext
    b = 2*transpose(Δuf)*(Δu + Δur) .+ 2*Δλ.*φ^2*transpose(fext)*fext
    c = transpose(Δu + Δur)*(Δu + Δur) .- Δl^2 .+ Δλ.^2.0.*φ^2.0*transpose(fext)*fext

    D = b.^2 .- (4).*a*c
    @assert D >= 0 "Discriminant must not be negative!"

    δλ = (-b .+ [-1;1]*sqrt(D))/2/a
    Δu1 = Δu + η*Δur + δλ[1]*Δuf
    Δu2 = Δu + η*Δur + δλ[2]*Δuf
    
    if transpose(Δu1)*Δu/Δl^2 >= transpose(Δu2)*Δu/Δl^2
        Δλ .= Δλ .+ δλ[1]
        Δu .= Δu1
    else
        Δλ .= Δλ .+ δλ[2]
        Δu .= Δu2
    end

    nothing
end

"""
# Ramm corrector step
[1] Fafard1993
"""
function rammcorrection!(Δu,Δλ,Kt,R,fext,Δl,(;Δur,Δuf,Δλbar);cylindrical=false)
    Δur .= -(Kt\R)
    Δuf .= Kt\fext

    Δλbar .= -(transpose(Δu)*Δur)/(transpose(Δu)*Δuf)
    
    Δu .+= Δur .+ Δλbar.*Δuf
    Δλ .+= Δλbar

    nothing
end

"""
# Modified Crisfield-Ramm corrector
[1] Fafard1993
"""
function mcrcorrection!(Δu,Δλ,Kt,R,fext,Δl,(;Δur,Δuf,Δλbar,α);cylindrical=false)
    rammcorrection!(Δu,Δλ,Kt,R,fext,Δl,(;Δur,Δuf,Δλbar);cylindrical=false)

    # correction to match the arc radius
    if cylindrical
        α .= sqrt(Δl^2)/(norm(Δu))
    else
        α .= sqrt(Δl^2)/(norm([Δu;Δλ]))
    end
    
    Δu .*= α
    Δλ .*= α

    nothing
end

end