
begin
    using QuadGK, Plots, NLsolve, CSV, DataFrames, ForwardDiff
    include("parameters.jl")
    include("functions.jl")
end

function gapsolver(mu,T,chuteinit)
    sistema = nlsolve(x->(dM(x[1],x[2],mu,T,x[3]),dphi(x[1],x[2],mu,T,x[3]),dphib(x[1],x[2],mu,T,x[3])),chuteinit)
    return sistema.zero
end

begin
    T_vals = range(1e-8, 0.35, length = 100)
    phi_vals = zeros(length(T_vals))
    mu = 0
    chuteinit = [0.01,0.01,0.3]
    Threads.@threads for i in 1:length(T_vals)
        T = T_vals[i]
        solution = gapsolver(mu,T)
        phi_vals[i] = solution[1]
        chuteinit = solution
    end
    plot(T_vals, phi_vals)
end

