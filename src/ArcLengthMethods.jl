module ArcLengthMethods

using NLsolve, LinearAlgebra, ForwardDiff

export arclengthmethod
export rikscorrection,crisfieldcorrection,rammcorrection,mcrcorrection

function arclengthmethod(fint,fext,correctorstep,Δl,u0;
    λ0=1e-2,ftol=1e-8,iterations=100,steps=1000,verbose=false,adaptivestep=true)
    
    
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
       
    step = 1
    while last(qs)[end] < 1
        
        u = last(qs)[1:end-1]
        λ = last(qs)[end]

        Δu,Δλ = secantpredictor(last(qs,2),Δl)

        Kt = ForwardDiff.jacobian(fint, u) 
        # Kt = lu(ForwardDiff.jacobian(fint, u)) #todo nefunguje u rikse
        
        iteration = 0
        R = similar(u).+1

        converged = false
        while iteration < iterations-1
            evalfun!(R,u + Δu,λ + Δλ)
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

    qs
end

function secantpredictor(qs,Δl)
    Δq = diff(qs)[1]
    Δq = Δl/norm(Δq)*Δq
    return Δq[1:end-1],last(Δq)
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
    δū = -(Kt\R)
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

Co takhle pouzit JFNK?
https://book.sciml.ai/notes/09/
znam smer ve kterém hledám derivaci...
"""
function rammcorrection(Δu,Δλ,Kt,R,fext,Δl)
    Δur = -(Kt\R)
    Δuf = Kt\fext

    Δλ2 = -(transpose(Δu)*Δur)/(transpose(Δu)*Δuf)
    
    Δu += Δur + Δλ2*Δuf
    Δλ += Δλ2 

    Δu,Δλ
end

"""
# Modified Crisfield-Ramm corrector
[1] Fafard1993
"""
function mcrcorrection(Δu,Δλ,Kt,R,fext,Δl)
    Δur = -(Kt\R)
    Δuf = Kt\fext

    Δλ2 = -(transpose(Δu)*Δur)/(transpose(Δu)*Δuf)
    
    Δu += Δur + Δλ2*Δuf
    Δλ += Δλ2 

    # correction to match the arc radius
    α = sqrt(Δl^2)/(norm([Δu;Δλ]))
    Δu *= α
    Δλ *= α

    Δu,Δλ
end

end