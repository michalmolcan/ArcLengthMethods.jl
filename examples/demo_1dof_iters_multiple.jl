using Plots,PlotlyJS, BenchmarkTools, Revise
plotlyjs()

includet("../src/ArcLengthMethods.jl")
using .ArcLengthMethods

function fint(a)
    θ₀ = pi/3
    return [(1/sqrt(1-2*a[1]*sin(θ₀) + a[1]^2) -1)*(sin(θ₀) - a[1])]
end

fext = [1.0]
Δl = 1e-1
u0 = [1e-6]

opts = (
    verbose = true,
    adaptivestep=false,
    cylindrical = false,
    save_iters = true
)

qs_riks,iters_data_riks = arclengthmethod(fint,fext,Δl,u0;method=:riks,opts...);
qs_crisfield,iters_data_crisfield = arclengthmethod(fint,fext,Δl,u0;method=:crisfield,opts...);
qs_ramm,iters_data_ramm = arclengthmethod(fint,fext,Δl,u0;method=:ramm,opts...);
qs_mcr,iters_data_mcr = arclengthmethod(fint,fext,Δl,u0;method=:mcr,opts...);

function plot_iters(qs_riks,iters_data_riks,label)

    iters_data_riks2 = mapreduce(x->[x[1],x[2][1]+x[4][1],x[3][1]+x[5][1],x[6][1],x[7][1]],hcat,iters_data_riks)'
    iters_data_riks3 = mapreduce(x->[x[1],x[2][1],x[3][1],x[4][1],x[5][1],x[6][1],x[7][1]],hcat,iters_data_riks)'
    converged_data_riks2 = reduce(hcat,qs_riks)
    plt = plot(size=(1200,750))
    # plt = plot()
    plot!(plt,converged_data_riks2[1,:],converged_data_riks2[2,:],label=false,xlabel="lambda",ylabel="u",title=label,color=:black)
    # savefig(plt,"riks_iters.html")

    new_data = []
    for (i,data_row) in enumerate(eachrow(iters_data_riks3))
        if i == 1
            push!(new_data,[data_row])
        elseif last(new_data[end])[1] >= data_row[1]
            push!(new_data,[data_row])
        else
            push!(new_data[end],data_row)
        end
    end

    function circleShape(h,k,r)
        θ = LinRange(0,2pi,1000)
        h .+ r*sin.(θ), k .+ r*cos.(θ)
    end

    for (i,data_row) in enumerate(new_data)
        (iters,u,λ,Δu,Δλ,Δl,R)=data_row[1]
        plot!(plt,circleShape(u,λ,Δl),linealpha=0.5,label=false,color=:orange)

        data3 = reduce(hcat,data_row)'
        plot!(plt,data3[:,2]+data3[:,4],data3[:,3]+data3[:,5],label=false,linestyle=:dash,color=:purple)

        plot!(plt,data3[1:1,2]+data3[1:1,4],data3[1:1,3]+data3[1:1,5],label=false,markershape=:cross,color=:red,seriestype=:scatter)
        if length(data_row) > 1
            plot!(plt,data3[2:end,2]+data3[2:end,4],data3[2:end,3]+data3[2:end,5],label=false,markershape=:cross,color=:purple,seriestype=:scatter)
        end

        if i > 1
            data_prev = new_data[i-1][end]
            p1 = [data_prev[2]+data_prev[4],data_prev[3]+data_prev[5]]
            p2 = [data3[1,2]+data3[1,4],data3[1,3]+data3[1,5]]
            pp = hcat(p1,p2)'

            plot!(plt,pp[:,1],pp[:,2],label=false,linestyle=:dash,linealpha=0.5,color=:red)
        end
    end
    plt
end

plt1 = plot_iters(qs_riks,iters_data_riks,"Riks")
plt2 = plot_iters(qs_crisfield,iters_data_crisfield,"Crisfield")
plt3 = plot_iters(qs_ramm,iters_data_ramm,"Ramm")
plt4 = plot_iters(qs_mcr,iters_data_mcr,"MCR")

savefig(plt1,"docs/files/riks_iters.html")
savefig(plt2,"docs/files/crisfield_iters.html")
savefig(plt3,"docs/files/ramm_iters.html")
savefig(plt4,"docs/files/mcr_iters.html")