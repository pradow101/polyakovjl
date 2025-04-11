
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
    T_vals = range(0.001, 0.3, length = 1000) #range for T
    phi_vals = zeros(length(T_vals)) # Arrays which will store the phi, phib and M solutions
    phib_vals = zeros(length(T_vals))
    M_vals = zeros(length(T_vals))
    mu = 0.34 
    chuteinit = [0.01,0.01,0.3]
     for i in 1:length(T_vals) #Initial guess
        T = T_vals[i]  #Tells the program to use the ith value of the T_vals array
        solution = gapsolver(mu, T, chuteinit) #Call gapsolver function, store it in the variable solution
        phi_vals[i] = solution[1] #solution is a vector of 3 floats, and we are storing the first one in phi_vals[i],
        phib_vals[i] = solution[2] #the second one in phib_vals[i], and the third one in M_vals[i]
        M_vals[i] = solution[3]
        chuteinit = solution #update the initial guess with the previous solution
    end
    plot(T_vals, [phi_vals, phib_vals], grid=true, gridalpha=0.5)
end
##debugging
begin
    dmvals = zeros(length(T_vals))
    dphivals = zeros(length(T_vals))
    for i in 1:length(T_vals)
        T = T_vals[i]
        dmvals[i] = abs(dM(phi_vals[i], phib_vals[i], mu, T_vals[i], M_vals[i]))
        dphivals[i] = abs(dphi(phi_vals[i], phib_vals[i], mu, T_vals[i], M_vals[i]))
    end
    df = DataFrame(T=T_vals,dM=dmvals)
    df2 = DataFrame(T=T_vals,dphi=dphivals)
    CSV.write("dphi.csv", df2)
    CSV.write("dM.csv", df)
    plot(T_vals, [dmvals,dphivals], label = ["dM" "dphi"], grid=true, gridalpha=0.5, xlabel = "T", ylabel = "dm/dphi")
end

begin
    Mi = range(0.01, 0.3, length = 100)
    phi_k = range(0.01, 0.3, length = 100)
    phib_k = range(0.01, 0.3, length = 100)
    yi = [potential(phi, 0.15, 0.34, 0.17, M) for (M, phi) in zip(Mi, phi_k)]
    plt = plot3d(Mi, phi_k, yi, xlabel = "M", ylabel = "phi", zlabel = "potential", title = "Potential vs M and phi")
end

begin
    # Define the range of M and phi values
    Mi = range(0.01, 0.3, length = 100)  # Range for M
    phi_k = range(0.01, 0.3, length = 100)  # Range for phi

    # Create a grid of M and phi values
    M_grid = repeat(Mi, 1, length(phi_k))  # Repeat M values along rows
    phi_grid = repeat(phi_k', length(Mi), 1)  # Repeat phi values along columns

    # Compute yi for each pair of (M, phi) on the grid
    yi = [potential(phi, 0.15, 0.34, 0.17, M) for (M, phi) in zip(M_grid[:], phi_grid[:])]

    # Reshape yi to match the grid dimensions
    yi = reshape(yi, length(Mi), length(phi_k))

    # Plot the 3D surface
    plot3d(Mi, phi_k, yi; xlabel = "M", ylabel = "phi", zlabel = "potential", title = "Potential vs M and phi")
end

begin
    # Define the range of M and phi values
    Mi = range(-0.5, 0.5, length = 100)  # Range for M
    phi_k = range(-0.5, 0.5, length = 100)  # Range for phi

    # Create a grid of M and phi values
    M_grid = repeat(Mi, 1, length(phi_k))  # Repeat M values along rows
    phi_grid = repeat(phi_k', length(Mi), 1)  # Repeat phi values along columns

    # Compute yi for each pair of (M, phi) on the grid
    yi = [potential(phi, 0.18, 0.34, 0.17, M) for (M, phi) in zip(M_grid[:], phi_grid[:])]

    # Reshape yi to match the grid dimensions
    yi = reshape(yi, length(Mi), length(phi_k))

    # Plot the 3D surface
    plot3d(Mi, phi_k, yi, st=:surface, xlabel = "M", ylabel = "phi", zlabel = "potential", title = "Potential vs M and phi")
end