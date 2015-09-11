function simWithUncer(theta,policyA1,EV,startA,e,logy1,randYpath)

# This function takes the policy functions and value functions, along with
# starting assets and returns simulated paths of income, consumption,
# assets and value

#% ------------------------------------------------------------------------
# Declare global we need this file have access to
#global mu sigma rho T r Tretire
#global Agrid  Ygrid numSims interpMethod uncertaintyMethod hcIncPDF;
#global normBnd

    paramsInco=theta[       thetaL[1]+1:sum(thetaL[1:2])];

  mu   =paramsInco[1];
  sigma=paramsInco[2];
  rho  =paramsInco[3];

#% ------------------------------------------------------------------------
# Initialise arrays that will hold the paths of income consumption, value
# and assets

# EV=exVal
# startingA=startA
# e=drawsNor
# logy1=logy1


# Arguments for output
y = Array(Float64,T, numSims);            # income
c = Array(Float64,T, numSims);            # consumption
v = Array(Float64,T, numSims);            # value
a = Array(Float64,T + 1,numSims);         # this is the path at the start of each period, so we include the 'start' of death

# Other arrays that will be used below
ly = Array(Float64,T, numSims);           # log income
ypathIndex = Array(Float64,T, numSims);   # holds the index (location) in the vector



      # *******************************************************************************************
      # Do all possible interpolations (get the objects) before any simulation, this will avoid the
      # repetition of this costly process.

      vinter   =Array(Any,T, numPointsY);
      ainter   =Array(Any,T, numPointsY);
      ixYf     =Array(Any,T);

      for t = 1:1:T

          A=Agrid[t]
          Y=vec(Ygrid[t,:])

          # Due to the odd restrictions of Grid package (even-spaced grid), we need this...
          incindi=1.:1:length(Ygrid[t,:])
          ixYf[t] = InterpIrregular(squeeze(Ygrid[t,:],1),[incindi],BCnearest, InterpNearest)

          for ixY = 1:1:numPointsY   # points on income grid

              # Given that, now we need to interpolate over A dimension
              tV =  squeeze(EV      [t, :,ixY],1);         #the relevant part of the value function
              tA1 = squeeze(policyA1[t, :,ixY],1);   #the relevant part of the policy function
              # --------------------------------------------------------------------
              vinter[t, ixY]  =CoordInterpGrid(A, tV,  BCreflect, InterpLinear)
              ainter[t, ixY]  =CoordInterpGrid(A, tA1, BCreflect, InterpLinear)
          end
      end

#% ------------------------------------------------------------------------
# Obtain paths using the initial condition and the policy and value
# functions
#-------------------------------------------------------------------------%#
# Scenario where we have hardcoded the discrete
# number of income draws and their pdf

if (uncertaintyMethod == 0)     # we have hardcoded the income process

    #First get random discrete draws using subroutine getDiscreteDraws.
    #More details given in that routine. Given a seed (which
    #we set equal to the time period), it draws randomly from vector Ygrid(t, :)
    #with pdf hcIncPDF
    for t = 1:1:T
        ypathIndex[t, :]=vec(randYpath[t,:])
        y[t, :] = Ygrid[t, vec(ypathIndex[t, :])];
    end

    for s = 1:1: numSims              # loop through individuals
        a[1, s] = startA;
        for t = 1:1:T                 # loop through time periods for a particular individual

            if (t >= Tretire)
                ixY = 1;
            else
                ixY = ypathIndex[t, s];
            end

              # Interpolate over policy rules
              ainterx  =ainter[t,ixY]
              vinterx  =vinter[t,ixY]

              a[t+1, s] = ainterx[ tLog(a[t, s],lbAGrid,gridMethod) ];
              v[t  , s] = vinterx[ tLog(a[t, s],lbAGrid,gridMethod) ];

            c[t, s] = a[t, s]  + y[t, s] - (a[t+1, s]/(1+r));
        end   #t
    end # s

#----------------------------------------%
# Scenario where income draws are normally distributed

 elseif (uncertaintyMethod == 1)

    # Initial conditions first
    sig_inc = sigma/ ((1-rho^2)^0.5);
    for s = 1:1: numSims
        a[1, s]  = startA;
        ly[1, s] = truncate(logy1[s], -normBnd*sig_inc , normBnd*sig_inc );
        y[1, s]  = exp(ly[1, s]);
    end

     # Now run their lifes!
     for t = 1:1:T                              # loop through time periods for a particular individual
      #t=40
         lbAGrid=borrowCon[t]; #Why it is independent of "j"? This is only costructed in order to be able to translate from A into the logA grid
         for s = 1:1: numSims

             if (t >= Tretire)                      # set income to zero if the individual has retired
                 y[t, s] = 0;
             end
            if (t < Tretire)                       # first for the before retirement periods
                #clear tA1 tV;                      # necessary as the dimensions of these change as we wor through this file

                tA1  = squeeze(policyA1[t, :, :],1);   # the relevant part of the policy function
                tV   = squeeze(EV[t, :, :],1);         # the relevant part of the value function


                #Get the income
                ixYfx=ixYf[t];
                ixY=ixYfx[y[t,s]];
                ypathIndex[t, s]=ixY;

             else                          # next for the post retirement periods

                #clear tA1 tV;
                tV  = EV[t, :, 1];         # the relevant part of the value function
                tA1 = policyA1[t, :, 1];  # the relevant part of the policy function

                ixY=1;

             end #% if (t < Tretire)

              # Interpolate over policy rules
              ainterx  =ainter[t,ixY]
              vinterx  =vinter[t,ixY]

              a[t+1, s] = ainterx[ tLog(a[t, s],lbAGrid,gridMethod) ];
              v[t  , s] = vinterx[ tLog(a[t, s],lbAGrid,gridMethod) ];


              if (t != T && t < (Tretire-1) )  # Get next year's income

                ly[t+1, s] = (1 -rho) * mu + rho * ly[t, s] + e[t + 1, s];
                ly[t+1, s] = truncate(ly[t+1, s], -normBnd*sig_inc,normBnd*sig_inc );
                y[t+1, s] = exp( ly[t+1, s] );
              elseif (t != T)
                y[t+1, s] = 0;

              end # if (t != T)


            # Check whether next period's asset is below the lowest
            # permissable
           # if ( a[t+1, s] < Agrid[t+1, 1] )
           #     a[t+1, s]  = checkSimExtrap( Agrid[t+1, 1],y[t, s], t );
           #  end
            c[t, s] = a[t, s]  + y[t, s] - (a[t+1, s]/(1+r));
        end   #s
     end # t
 end # elseif (uncertaintyMethod == 1)
# %-------------------------------------------------------------------------%
# %-------------------------------------------------------------------------%


 return  y, c, a, v

end
