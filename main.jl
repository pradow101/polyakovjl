
begin
    using QuadGK, Plots, NLsolve, CSV, DataFrames, ForwardDiff
    include("parameters.jl")
    include("functions.jl")
end

function gapsolver(mu,T)
    sistema = nlsolve(x->(dM(x[1],x[2],mu,T,x[3]),dphi(x[1],x[2],mu,T,x[3]),dphib(x[1],x[2],mu,T,x[3])), [0.01,0.01,0.3]).zero
end


begin
    M_vals = range(0,0.5,length=50)
    yi = [dM(1e-8, 1e-8, 0.2, 0.01, Mi) for Mi in M_vals]
    plot(M_vals, yi, label="dM", xlabel="M", ylabel="dM", title="dM vs M", grid=true, gridalpha=0.5)
end
