begin
    using QuadGK, Plots, NLsolve, CSV, DataFrames, ForwardDiff, DataInterpolations, LocalFunctionApproximation, Interpolations
    include("parameters.jl")
    include("functions.jl")
    include("main.jl")

    plotly()
end

##Potencial no crossover
begin
    ## Ω x M
    Mvcrossover = range(-0.6,0.6,100)
    Pvcrossover = [potential(0.1,0.1,0.2,0.2, M) for M in Mvcrossover]
    plot(Mvcrossover, Pvcrossover)
end

##Potencial no ponto crítico
begin
    ## Ω x M
    Mvcritical = range(-0.6,0.6,100)
    Pvcritical = [potential(0.1,0.1,0.2,0.2, M) for M in Mvcritical]
    plot(Mvcritical, Pvcritical)
end

##Potencial na primeira ordem
begin
    ## Ω x M
    Mvfirst = range(-0.6,0.6,100)
    Pvfirst = [potential(0.1,0.1,0.2,0.2, M) for M in Mvfirst]
    plot(Mvfirst, Pvfirst)
end

