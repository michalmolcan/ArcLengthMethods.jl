using Plots, BenchmarkTools, Revise

includet("../src/ArcLengthMethods.jl")
using .ArcLengthMethods

function fint(a)
    θ₀ = pi/3
    return [(1/sqrt(1-2*a[1]*sin(θ₀) + a[1]^2) -1)*(sin(θ₀) - a[1])]
end


function demo()
    fext = [1.0]
    Δl = 5e-2
    u0 = [1e-6]

    verbose = false

    # Riks arc length method
    println("Riks arc length method")
    qs = arclengthmethod(fint,fext,Δl,u0;verbose=verbose,method=:riks)
    println("Riks method done")
    pl = plot([u[1] for u in qs],[u[2] for u in qs],ls=:auto,legend=:bottomright,label="Riks")
    ylabel!("load factor λ")
    xlabel!("displacement u₂")

    # # Crisfields arc length method
    # println("Crisfields arc length method")
    # qs = arclengthmethod(fint,fext,Δl,u0,verbose=verbose,method=:crisfield)
    # println("Crisfields method done")
    # plot!([u[1] for u in qs],[u[2] for u in qs],ls=:auto,label="Crisfield")

    # # Ramm arc length method
    # println("Ramms arc length method")
    # qs = arclengthmethod(fint,fext,Δl,u0;verbose=verbose,method=:ramm)
    # println("Ramms method done")
    # plot!([u[1] for u in qs],[u[2] for u in qs],ls=:auto,label="Ramm")

    # # MCR arc length method
    # println("MCR arc length method")
    # qs = arclengthmethod(fint,fext,Δl,u0;verbose=verbose,method=:mcr)
    # println("MCR method done")
    # plot!([u[1] for u in qs],[u[2] for u in qs],ls=:auto,label="MCR")

    return pl
end

demo()