using DifferentialEquations
using Random
using NPZ
#============================================================#
function int_to_trinary(a)
    s = "" ;
    quotient = a ;
    while quotient > 2
        rem = quotient % 3 ;
        quotient = div(quotient,3) ;
        s = string(rem) * s ;
    end
    s = string(quotient) * s ;
    if length(s) < 6
        for i = 1:(6-length(s))
            s = '0'*s ;
        end
    end
    return(s) ;
end
#=============================================================#
function get_omegas(ctstring,theta4,theta5,theta6)
    if ctstring[4] == '0'
        omegaA = 0 ;
        omegaD = 0 ;
    elseif ctstring[4] == '1'
        omegaA = theta4 ;
        omegaD = 1 ;
    elseif ctstring[4] == '2'
        omegaA = 0 ;
        omegaD = theta4 ;
    end
    #-----------------------#
    if ctstring[5] == '0'
        omegaB = 0 ;
        omegaE = 0 ;
    elseif ctstring[5] == '1'
        omegaB = theta5 ;
        omegaE = 1 ;
    elseif ctstring[5] == '2'
        omegaB = 0 ;
        omegaE = theta5 ;
    end
    #-----------------------#
    if ctstring[6] == '0'
        omegaC = 0 ;
        omegaF = 0 ;
    elseif ctstring[6] == '1'
        omegaC = theta6 ;
        omegaF = 1 ;
    elseif ctstring[6] == '2'
        omegaC = 0 ;
        omegaF = theta6 ;
    end
    #-------------------------#
    return(omegaA,omegaB,omegaC,omegaD,omegaE,omegaF) ;
end

#=============================================================#
function get_final_state(sol,NCells)
    A = sol[(1*NCells+1):2*NCells,end] ;
    B = sol[(2*NCells+1):3*NCells,end] ;
    C = sol[(3*NCells+1):4*NCells,end] ;
    States = zeros(Int8,NCells) ;
    for i = 1:NCells
        if (A[i] > B[i]) & (A[i] > C[i])
            States[i] = 0 ;
        elseif (B[i] > A[i]) & (B[i] > C[i])
            States[i] = 1 ;
        else
            States[i] = 2 ;
        end
    end
    return(States)
end
#=====================================================#
function coupled_frenchflag_ultimate!(du,u,p,t)
    Kappa,h1,h2,h3,h4,h5,k1,k2,k3,alpha,beta,gamma,tauD,tauNb,betaD,betaNb,K_N,NCells,Smax,Slambda,theta1,theta2,theta3,theta4,theta5,theta6,ctstring = p ;

    S = u[1:NCells]
    A  = u[(1*NCells+1):2*NCells]
    B  = u[(2*NCells+1):3*NCells]
    C  = u[(3*NCells+1):4*NCells]
    D  = u[(4*NCells+1):5*NCells]
    Nb = u[(5*NCells+1):6*NCells]
    for i in 1:NCells
        if i == 1
            D_N = D[2] ;
        elseif i == NCells
            D_N = D[NCells-1] ;
        else
            D_N = 0.5 * (D[i-1] + D[i+1]) ;
        end
        #------------------------------------#
        du[i] = (Smax* exp(-Slambda*(i-1))) - S[i] ;
        du[(5*NCells)+i] = betaNb * (D_N/(1 + D_N)) - Nb[i]/tauNb ;
        #-------------------------------------#
        if ctstring[1] == '0'
            du[(1*NCells)+i] = (alpha)/(1 + (C[i]/Kappa)^h1 + (B[i]/Kappa)^h2) - k1*A[i] ;
        elseif ctstring[1] == '1'
            du[(1*NCells)+i] = (alpha + theta1*Nb[i])/(1 + (C[i]/Kappa)^h1 + (B[i]/Kappa)^h2 + Nb[i]) - k1*A[i] ;
        elseif ctstring[1] == '2'
            du[(1*NCells)+i] = (alpha)/(1 + (C[i]/Kappa)^h1 + (B[i]/Kappa)^h2 + theta1*Nb[i]) - k1*A[i] ;
        end
        #--------------------------------------#
        if ctstring[2] == '0'
            du[(2*NCells)+i] = ((beta*S[i])/(1 + S[i]))*(1/(1 + (C[i]/Kappa)^h3))  - k2*B[i] ;
        elseif ctstring[2] == '1'
            du[(2*NCells)+i] = ((beta*S[i] + theta2*Nb[i])/(1 + S[i] + Nb[i]))*(1/(1 + (C[i]/Kappa)^h3))  - k2*B[i] ;
        elseif ctstring[2] == '2'
            du[(2*NCells)+i] = ((beta*S[i])/(1 + S[i] + theta2*Nb[i]))*(1/(1 + (C[i]/Kappa)^h3))  - k2*B[i] ;
        end
        #--------------------------------------#
        if ctstring[3] == '0'
            du[(3*NCells)+i] = ((gamma*S[i])/(1 + S[i]))*(1/(1 + (B[i]/Kappa)^h4 + (A[i]/Kappa)^h5)) - k3*C[i] ;
        elseif ctstring[3] == '1'
            du[(3*NCells)+i] = ((gamma*S[i] + theta3*Nb[i])/(1 + S[i] + Nb[i]))*(1/(1 + (B[i]/Kappa)^h4 + (A[i]/Kappa)^h5)) - k3*C[i] ;
        elseif ctstring[3] == '2'
            du[(3*NCells)+i] = ((gamma*S[i])/(1 + S[i] + theta3*Nb[i]))*(1/(1 + (B[i]/Kappa)^h4 + (A[i]/Kappa)^h5)) - k3*C[i] ;
        end
        #---------- remaining three strings to be considered together giving raise to 8 possibilities ----------------------------#
        omegaA,omegaB,omegaC,omegaD,omegaE,omegaF = get_omegas(ctstring,theta4,theta5,theta6) ;
        du[(4*NCells)+i] = (betaD + omegaA*A[i] + omegaB*B[i] + omegaC*C[i] )/(1 + omegaD*A[i] + omegaE*B[i] + omegaF*C[i] ) - D[i]/tauD ;

    end
end#function
#==================================================================================================#
#CID = 1 : 728
function Explore_All_Thetas(CID,niters)
    PATH = "Data/" ;

    n_vals = 100 ;

    ctstring = int_to_trinary(CID) ;
    println(ctstring) ;

    NCells = 30 ;
    Smax,Slambda = 100, 0.3 ;
    Kappa = 1.0 ;
    h1,h2,h3,h4,h5 = 6,2,5,1,1 ;
    k1,k2,k3 = 1.0,1.0,1.0 ;
    alpha,beta,gamma = 4.0,6.3,5.0 ;
    tauD,tauNb = 1.0,1.0 ;
    betaD,betaNb = 5.0,5.0 ;
    K_N = 1.0 ;
    OD = zeros(Int16,niters,NCells) ;

    Thetas = zeros(6) ;
    ThetaUP = LinRange(1.0,10.0,n_vals);
    ThetaDN = LinRange(0.1,1.0,n_vals);

    for i = 1:niters
        println(i) ;
        for j = 1:6
            if ctstring[j] == '0'
                Thetas[j] = 0.0 ;
            elseif ctstring[j] == '1'
                Thetas[j] = ThetaUP[rand(1:n_vals)] ;
            elseif ctstring[j] == '2'
                Thetas[j] = ThetaDN[rand(1:n_vals)] ;
            end
        end
        theta1,theta2,theta3,theta4,theta5,theta6 = Thetas ;
        #--------------------------------------------------#

        p = Kappa,h1,h2,h3,h4,h5,k1,k2,k3,alpha,beta,gamma,tauD,tauNb,betaD,betaNb,K_N,NCells,Smax,Slambda,theta1,theta2,theta3,theta4,theta5,theta6,ctstring ;
        u0 = zeros(6*NCells) ;
        for k = 1 : NCells
            u0[k] = (Smax* exp(-Slambda*(k-1))) ;
        end

        tmax = 100 ;
        tspan = (0.0,tmax) ;
        prob  = ODEProblem(coupled_frenchflag_ultimate!,u0,tspan,p) ;
        sol   = solve(prob) ;
        FS = get_final_state(sol,NCells) ;
        OD[i,:] = FS ;
    end
    ofilename = PATH * "Ensemble_Ultimate_NCells_"*string(NCells)*"_CID_"*string(CID)*".npy" ;
    npzwrite(ofilename,OD) ;
end#function
#===================================================================================#


CID = 345 ;
niters = 100 ;
Explore_All_Thetas(CID,niters) ;
