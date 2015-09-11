function solveValueFunction(theta::Array{Float64,1})

#This function obtains the value function for each time period and
#the policy function (i.e. optimal next-period asset choice) for each time
#period. From there we can work the optimal consumption level.

#The approach taken is by backwards recursion. The optimisation each period
#is carried out using 'fminbnd'. This is an in-built optimiser in Matlab.
#The optimisation routine it uses is known as the 'golden search method'

  paramsUtil=theta[                 1:thetaL[1]];
  paramsInco=theta[       thetaL[1]+1:sum(thetaL[1:2])];

  mu   =paramsInco[1];
  sigma=paramsInco[2];
  rho  =paramsInco[3];

#% ------------------------------------------------------------------------
# GENERATE MATRICES TO STORE NUMERICAL APPROXIMATIONS AND INITIATE AS NAN

# Matrices to hold the policy, value and marginal utility functions
V        = Array(Float64,T+1, numPointsA, numPointsY);
policyA1 = Array(Float64,T,   numPointsA, numPointsY);
policyC  = Array(Float64,T,   numPointsA, numPointsY);

#Matrices to hold expected value and marginal utility functions
EV  = Array(Float64,T+1, numPointsA, numPointsY);


#% ------------------------------------------------------------------------
#Set the terminal value function and expected value function to 0

EV[T + 1, :,:]  = .0;          # continuation value at T-1
V[T + 1,:,:]    = .0;

lbAGrid=borrowCon[T+1];

#% ------------------------------------------------------------------------
# SOLVE RECURSIVELY THE CONSUMER'S PROBLEM, STARTING AT TIME T-1 AND MOVING
# BACKWARDS TO ZERO, ONE PERIOD AT A TIME


for ixt=T:-1:1                                # Loop from time T-1 to 1
#ixt=T

    Agrid1=Agrid[ixt + 1]

    lbAGrid=borrowCon[ixt];
    lbA1 = eExp( Agrid[ixt + 1][1] , lbAGrid,gridMethod);          # lower bound: assets tomorrow

    Agrido=[Agrid[ixt]]     # This period A grid
    Ygrido=Ygrid[ixt, :]    # This period Y grid

    for ixA = 1:1:numPointsA                  # points on asset grid
#ixA=1

        # STEP 1. solve problem at grid points in assets and income
        # ---------------------------------------------------------
        A  = eExp( Agrido[ixA] ,lbAGrid,gridMethod);            # assets today


        for ixY = 1:1:numPointsY               # points on income grid

        #ixY=1

            # Value of income and information for optimisation
            Y    = Ygrido[ ixY];            # income today
            ubA1 = (A + Y - minCons)*(1+r);    # upper bound: assets tomorrow
            EV1  = vec(EV[ixt + 1,:, ixY]);    # relevant section of EV matrix (in assets tomorrow), vector

            intfV  = CoordInterpGrid(Agrid1,EV1, BCreflect, InterpLinear); # Define the interpolation

            # Compute solution
            if (ubA1 - lbA1 < minCons)                       # if liquidity constrained
                negV = objectivefunc(lbA1,[0],paramsUtil,r,minCons,gridMethod,lbAGrid,A, Y,intfV);
                policyA1[ixt,ixA,ixY] = lbA1;
            else                                             # if interior solution

                # Using NLopt...
                opt = Opt(:LN_SBPLX, 1)
                lower_bounds!(opt, [lbA1])
                upper_bounds!(opt, [ubA1])
                xtol_rel!(opt,1e-4)
                min_objective!(opt, (a1,grad::Vector) -> objectivefunc(a1,grad,paramsUtil,r,minCons,gridMethod,lbAGrid, A, Y,intfV))
                (valor, ax, toler)=optimize(opt,[ ((ubA1-abs(lbA1))/2) ])

                policyA1[ixt,ixA,ixY]=ax[1];
                negV=valor;

                if isnan(negV)
                    println("upps, there is an error with the optimization in ixt=$ixt, ixA=$ixA, ixY=$ixY")
                end

            end # if (ubA1 - lbA1 < minCons)

            # Store solution and its value
            policyC[ixt, ixA, ixY] = A + Y - policyA1[ixt, ixA, ixY]/(1+r);
            V[ixt, ixA, ixY]       = -negV;
        end #ixY


        # STEP 2. integrate out income today conditional on income
        # yesterday to get EV and EdU
        # --------------------------------------------------------
        realisedV = squeeze(V[ixt, ixA, :],1);
        for ixY = 1:1:numPointsY
            d=incTransitionMrx[ixY,:]*realisedV'
            EV[ixt, ixA, ixY]  = d[1,1];
        end #ixY
    end #ixA

    println("Passed period $ixt of $T")

end #ixt

return policyA1, policyC, V, EV

end #function

