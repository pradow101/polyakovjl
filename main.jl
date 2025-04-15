
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
    T_vals = range(0.1, 0.4, length = 100) #range for T
    phi_vals = zeros(length(T_vals)) # Arrays which will store the phi, phib and M solutions
    phib_vals = zeros(length(T_vals))
    M_vals = zeros(length(T_vals))
    mu = 0
    chuteinit = [0.01,0.01,0.3]
     for i in 1:length(T_vals) #Initial guess
        T = T_vals[i]  #Tells the program to use the ith value of the T_vals array
        solution = gapsolver(mu, T, chuteinit) #Call gapsolver function, store it in the variable solution
        phi_vals[i] = solution[1] #solution is a vector of 3 floats, and we are storing the first one in phi_vals[i],
        phib_vals[i] = solution[2] #the second one in phib_vals[i], and the third one in M_vals[i]
        M_vals[i] = solution[3]
        chuteinit = solution #update the initial guess with the previous solution
    end
    plot(T_vals, [M_vals, phi_vals,phib_vals], grid=true, gridalpha=0.5, xlabel = "T", ylabel = "phi, M", title = "M and phi solutions")
end

        
begin       #calculating and plotting the pressure
    pf_vals = zeros(length(T_vals))
    for i in 1:length(T_vals)
        T = T_vals[i]
        phi = phi_vals[i]
        phib = phib_vals[i]
        M = M_vals[i]
        pf_vals[i] = -(potential(phi, phib, 0, T, M) - potential(phi, phib, 0, 0.001, M))/pf(T) #pressure
    end
    plot(T_vals, pf_vals, grid = true, gridalpha=0.5, xlabel = "T", ylabel = "Pressure", title = "Pressure vs T")
end

begin
    
end



##debugging
# begin
#     dmvals = zeros(length(T_vals))
#     dphivals = zeros(length(T_vals))
#     for i in 1:length(T_vals)
#         T = T_vals[i]
#         dmvals[i] = abs(dM(phi_vals[i], phib_vals[i], mu, T_vals[i], M_vals[i]))
#         dphivals[i] = abs(dphi(phi_vals[i], phib_vals[i], mu, T_vals[i], M_vals[i]))
#     end
#     df = DataFrame(T=T_vals,dM=dmvals)
#     df2 = DataFrame(T=T_vals,dphi=dphivals)
#     CSV.write("dphi.csv", df2)
#     CSV.write("dM.csv", df)
#     plot(T_vals, [dmvals,dphivals], label = ["dM" "dphi"], grid=true, gridalpha=0.5, xlabel = "T", ylabel = "dm/dphi")
# end


