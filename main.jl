begin
    using QuadGK, Plots, NLsolve, CSV, DataFrames, ForwardDiff, Interpolations, LocalFunctionApproximation, CurveFit

    include("parameters.jl")
    include("functions.jl")

    plotly()
end

function gapsolver(mu,T,chutealto)
    sistema1 = nlsolve(x->(dM(x[1],x[2],mu,T,x[3]),dphi(x[1],x[2],mu,T,x[3]),dphib(x[1],x[2],mu,T,x[3])),chutealto)
    return sistema1.zero
end


function maxfind(y, T) #Here x will play the role of the derivatives of phi, phib and M solutions for a given μ
    for i in 150:length(y)-1
        if y[i+1] < y[i] && y[i-1] < y[i]  
            return T[i], y[i]
        end
    end
    return NaN, NaN # Return NaN if no maximum is found
end

begin
    function Trangesolver(mu, T_vals)
        phi_vals = zeros(length(T_vals)) # Arrays which will store the phi, phib and M solutions
        phib_vals = zeros(length(T_vals))
        M_vals = zeros(length(T_vals))
        chutealto = [0.01,0.01,0.4]
        for i in 1:length(T_vals) #Initial guess
            T = T_vals[i]  #Tells the program to use the ith value of the T_vals array
            solution = gapsolver(mu, T, chutealto) #Call gapsolver function, store it in the variable solution
            phi_vals[i] = solution[1] #solution is a vector of 3 floats, and we are storing the first one in phi_vals[i],
            phib_vals[i] = solution[2] #the second one in phib_vals[i], and the third one in M_vals[i]
            M_vals[i] = solution[3]
            #here I need to implement a condition to switch initial guess values
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
            derinterp[i,2] = only(Interpolations.gradient(itpphi, T_vals[i]))
            derinterp[i,3] = only(Interpolations.gradient(itpphib, T_vals[i]))
            derinterp[i,4] = -only(Interpolations.gradient(itpM, T_vals[i]))  
        end
        return derinterp
    end

    function murangesolver(T_vals)
        mu_vals = range(0,0.32,length=100)
        solutions = zeros(length(mu_vals), length(T_vals), 4)
        println(size(solutions))
        println(size(mu_vals))

        Threads.@threads for i in eachindex(mu_vals)
            solsi = Trangesolver(mu_vals[i], T_vals)
            #println(size(solsi))
            solutions[i,:,1] = solsi[1]
            solutions[i,:,2] = solsi[2]
            solutions[i,:,3] = solsi[3]
            solutions[i,:,4] = solsi[4]
        end
        return solutions, mu_vals
    end
end


# begin
#     T_vals = range(0.01,0.1,150)
#     mu = 0.34
#     T_vals, phi_vals, phib_vals, M_vals = Trangesolver(mu, T_vals)
#     plot(T_vals, [M_vals,phi_vals, phib_vals], grid=true, gridalpha=0.5, xlabel = "T", ylabel = "phi, phib", title = "phi and phib solutions", linewidth = 2, label = ["M" "ϕ" "ϕ*"])

# end

@time begin
    T_vals = range(0.04,0.4,1000)
    murange, muvalores = murangesolver(T_vals)
end

begin #ok, estão certos os resultados. Preciso achar um jeito de interpolar cada um
    a = murange[1,:,1]
    b = murange[1,:,2]
    c = murange[50,:,2]
    d = murange[100,:,2]
    e = murange[1,:,4]*3
    f = murange[50,:,4]*3
    g = murange[100,:,4]*3

    df = DataFrame(a = a, b = b, c = c, d = d, e = e, f = f, g = g)
    CSV.write("output.csv", df)
end

begin
    println(muvalores)
end

begin
    T_valores = zeros(length(murange[1,:,1]),length(muvalores))
    phi_valores = zeros(length(murange[1,:,1]),length(muvalores))
    phib_valores = zeros(length(murange[1,:,1]),length(muvalores))
    M_valores = zeros(length(murange[1,:,1]),length(muvalores))
    for i in eachindex(muvalores)
        T_vals, phi_vals, phib_vals, M_vals = murange[i,:,1], murange[i,:,2], murange[i,:,3], murange[i,:,4]
        interploop = Interp(T_vals, phi_vals, phib_vals, M_vals)
        T_valores[:,i] = interploop[:,1]
        phi_valores[:,i] = interploop[:,2]
        phib_valores[:,i] = interploop[:,3]
        M_valores[:,i] = interploop[:,4]
    end
end 



#plot(T_vals, [M_vals], grid=true, gridalpha=0.5, xlabel = "T", ylabel = "phi, M", title = "M and phi solutions")end

        
# begin       #calculating and plotting the pressure
#     T_vals = range(0.1,0.4,500)
#     mu = 0
#     T_vals, phi_vals, phib_vals, M_vals = Trangesolver(mu, T_vals)
#     pf_vals = zeros(length(T_vals))
#     for i in 1:length(T_vals)
#         T = T_vals[i]
#         phi = phi_vals[i]
#         phib = phib_vals[i]
#         M = M_vals[i]
#         pf_vals[i] = -(potential(phi, phib, 0, T, M) - potential(phi, phib, 0, 0.001, M))/pf(T) #pressure
#     end
#     plot(T_vals, pf_vals, grid = true, gridalpha=0.5, xlabel = "T", ylabel = "Pressure", title = "Pressure vs T", xrange = (0.1,0.395))
# end


begin
    Ttransitionphi = zeros(length(muvalores))
    Ttransitionphib = zeros(length(muvalores))
    actualphi = zeros(length(muvalores))
    TtransitionM = zeros(length(muvalores))
    Mutransition = muvalores
    Threads.@threads for i in 1:length(muvalores)
        Ttransitionphi[i] = maxfind(phi_valores[:,i], T_valores[:,1])[1]
        Ttransitionphib[i] = maxfind(phib_valores[:,i], T_valores[:,1])[1]
        TtransitionM[i] = maxfind(M_valores[:,i], T_valores[:,1])[1]
        actualphi[i] = (Ttransitionphi[i] + Ttransitionphib[i])/2
    end
    plot(Mutransition, [Ttransitionphi, actualphi], 
    label = ["ϕ Transition" "M Transition"],
    xlabel = "μ [GeV]",
    ylabel = "T [GeV]",
    title = "PNJL Phase Diagram", dpi=800, linewidth = 3)
    plot!([0.33179511923942306],[0.06557512531076132],color="green")
end


##Aqui vou definir outro jeito de se obter os valores críticos de μ usando a densidade
function gapsolvedensidade(T, chuteinit, nb)
    sistema = nlsolve(x->(dM(x[1],x[2],x[3],T,x[4]),dphi(x[1],x[2],x[3],T,x[4]),dphib(x[1],x[2],x[3],T,x[4]),densidade(x[1],x[2],x[3],T,x[4],nb)),chuteinit)
    return sistema.zero
end


##ESSA SEÇÃO DE CÓDIGO É A MAIS IMPORTANTE ATÉ AGORA, NÃO QUEBRAR
begin
    Nbvals = range(0.0001,0.01,length=1000)
    T = 0.03
    phi_vals = zeros(length(Nbvals)) # Arrays which will store the phi, phib and M solutions
    phib_vals = zeros(length(Nbvals))
    M_vals = zeros(length(Nbvals))
    mu_vals = zeros(length(Nbvals))
    potential_vals = zeros(length(Nbvals))
    chuteinit = [0.01,0.01,0.4,0.4]
    for i in 1:length(Nbvals) #Initial guess
        nb = Nbvals[i]  #Tells the program to use the ith value of the T_vals array
        solution = gapsolvedensidade(T, chuteinit, nb) #Call gapsolver function, store it in the variable solution
        phi_vals[i] = solution[1] #solution is a vector of 3 floats, and we are storing the first one in phi_vals[i],
        phib_vals[i] = solution[2] #the second one in phib_vals[i], and the third one in M_vals[i]
        mu_vals[i] = solution[3]
        M_vals[i] = solution[4]
        chuteinit = solution
        potential_vals[i] = potential(phi_vals[i], phib_vals[i], mu_vals[i], T, M_vals[i])
    end
    scatter(mu_vals, [M_vals], seriestype=:path)
end

begin
    #aqui, preciso dar um jeito de achar o valor para qual o potencial começa a voltar
    #vou usar a condição de quando o valor de mu aumenta pela primeira vez, de depois quando volta a cair.
    firstcurve = zeros(length(mu_vals))
    secondcurve = zeros(length(mu_vals))
    i=2
    while abs(potential_vals[i-1]) > abs(potential_vals[i])
        firstcurve[i] = potential_vals[i]
        i += 1
    end
    println(firstcurve)
end


begin
    interp1potenvals = interpolate((mu_vals,), potential_vals, Gridded(Linear()))
    interpolado = zeros(length(mu_vals))
    for i in range(length(mu_vals))
        interpolado[i] = only(Interpolations.gradient(interp1potenvals, mu_vals[i]))
    end
    plot(interpolado[i], mu_vals)
end