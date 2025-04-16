begin
    using QuadGK, Plots, NLsolve, CSV, DataFrames, ForwardDiff, Interpolations, LocalFunctionApproximation
    include("parameters.jl")
    include("functions.jl")

    plotly()
end

function gapsolver(mu,T,chutealto)
    sistema1 = nlsolve(x->(dM(x[1],x[2],mu,T,x[3]),dphi(x[1],x[2],mu,T,x[3]),dphib(x[1],x[2],mu,T,x[3])),chutealto)
    return sistema1.zero
end


function maxfind(x, T_vals) #Here x will play the role of the derivatives of phi, phib and M solutions for a given μ
    for i in 1:length(x)
        if x[i+1] < x[i] && x[i-1] < x[i]  
            return x[i], T_vals[i]
        end
    end
end

begin
    function Trangesolver(mu, T_vals)
        T_vals = range(0.01, 0.35, length = 500) #range for T
        phi_vals = zeros(length(T_vals)) # Arrays which will store the phi, phib and M solutions
        phib_vals = zeros(length(T_vals))
        M_vals = zeros(length(T_vals))
        mu = 0.27
        chutealto = [0.01,0.01,0.4]
        for i in 1:length(T_vals) #Initial guess
            T = T_vals[i]  #Tells the program to use the ith value of the T_vals array
            solution = gapsolver(mu, T, chutealto) #Call gapsolver function, store it in the variable solution
            phi_vals[i] = solution[1] #solution is a vector of 3 floats, and we are storing the first one in phi_vals[i],
            phib_vals[i] = solution[2] #the second one in phib_vals[i], and the third one in M_vals[i]
            M_vals[i] = solution[3]
            chutealto = solution
        end
        return T_vals, phi_vals, phib_vals, M_vals
    end

    function Interp(T_vals, phi_vals, phib_vals, M_vals)
        itpM = interpolate((T_vals,), M_vals, Gridded(Linear()))
        itpphi = interpolate((T_vals,), phi_vals, Gridded(Linear()))
        itpphib = interpolate((T_vals,), phib_vals, Gridded(Linear()))
        interp = zeros(length(T_vals), 4)
        derinterp = zeros(length(T_vals), 4)
        interp[:,1] = T_vals
        derinterp[:,1] = T_vals
        for i in range(1, length(T_vals))
            interp[i,2] = itpM(T_vals[i])
            interp[i,3] = itpphi(T_vals[i])
            interp[i,4] = itpphib(T_vals[i])
            derinterp[i,2] = only(Interpolations.gradient(itpphi, T_vals[i]))/23
            derinterp[i,3] = only(Interpolations.gradient(itpphib, T_vals[i]))/23
            derinterp[i,4] = -only(Interpolations.gradient(itpM, T_vals[i]))/23    
        end
        return interp, derinterp
    end 
    #interp, derinterp = Interp(T_vals, phi_vals, phib_vals, M_vals)
end




#plot(T_vals, [M_vals], grid=true, gridalpha=0.5, xlabel = "T", ylabel = "phi, M", title = "M and phi solutions")end

        
begin       #calculating and plotting the pressure
    pf_vals = zeros(length(T_vals))
    for i in 1:length(T_vals)
        T = T_vals[i]
        phi = phi_vals[i]
        phib = phib_vals[i]
        M = M_vals[i]
        pf_vals[i] = -(potential(phi, phib, 0, T, M) - potential(phi, phib, 0, 0.001, M))/pf(T) #pressure
    end
    plot(T_vals, pf_vals, grid = true, gridalpha=0.5, xlabel = "T", ylabel = "Pressure", title = "Pressure vs T", xrange = (0.1,0.395))
end



