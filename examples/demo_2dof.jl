using Plots, BenchmarkTools, Revise

includet("../src/ArcLengthMethods.jl")
using .ArcLengthMethods


"""
Example of 2dof beam construction from Vasios2015
"""
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

function demo1()

    # problem definition
    fext = [0, 1] 
    Δl = 5e-2
    u0 = [0,1e-6]

    # Riks arc length method
    println("Riks arc length method")
    qs = arclengthmethod(fint,fext,Δl,u0;method=:riks,verbose=true,adaptivestep=false)
    println("Riks method done")
    pl = plot([u[2] for u in qs],[u[3] for u in qs],ls=:auto,legend=:bottomright,label="Riks")
    ylabel!("load factor λ")
    xlabel!("displacement u₂")

    # Crisfields arc length method
    println("Crisfields arc length method")
    qs = arclengthmethod(fint,fext,Δl,u0;method=:crisfield,verbose=true,adaptivestep=false)
    println("Crisfields method done")
    plot!([u[2] for u in qs],[u[3] for u in qs],ls=:auto,label="Crisfield")

    # Ramm arc length method
    println("Ramms arc length method")
    qs = arclengthmethod(fint,fext,Δl,u0;method=:ramm,verbose=true,adaptivestep=false)
    println("Ramms method done")
    plot!([u[2] for u in qs],[u[3] for u in qs],ls=:auto,label="Ramm")

    # MCR arc length method
    println("MCR arc length method")
    qs = arclengthmethod(fint,fext,Δl,u0;method=:mcr,verbose=true,adaptivestep=false)
    println("MCR method done")
    plot!([u[2] for u in qs],[u[3] for u in qs],ls=:auto,label="MCR")

    return pl
end

demo1()

function demobenchmark()
    fext = [0, 1] 
    Δl = 5e-2
    u0 = [0,1e-6]
    # Benchmark 
    @btime arclengthmethod($fint,$fext,$Δl,$u0;method=:riks,adaptivestep=false)
    @btime arclengthmethod($fint,$fext,$Δl,$u0;method=:crisfield,adaptivestep=false)
    @btime arclengthmethod($fint,$fext,$Δl,$u0;method=:ramm,adaptivestep=false)
    @btime arclengthmethod($fint,$fext,$Δl,$u0;method=:mcr,adaptivestep=false)

    nothing
end

demobenchmark()