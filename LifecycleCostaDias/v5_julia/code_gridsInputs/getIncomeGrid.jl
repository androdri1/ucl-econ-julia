function getIncomeGrid(paramsInco)

# A function that returns:
# 1. an income grid
# 2. A Markovian transition matrix (Q) over income realisations
# 3. A vector of minimum incomes in each year
# 4. A vector of maximum incomes in each year

# ------------------------------------------------------------------------
#  Declare the global variables that will be necessary in this function
#global T Tretire sigma mu rho isUncertainty        % Structural economic paramters

mu    = paramsInco[1];
sigma = paramsInco[2];
rho   = paramsInco[3];

#global hcIncome hcIncPDF                        % Income process (hardcoded)
#global normBnd numPointsY uncertaintyMethod     % Numerical methods

# ------------------------------------------------------------------------


#----------------------------------------%
# Scenario where there is uncertainty - income draws are log normally distributed
#----------------------------------------%

#First get the standard deviation of income (from sigma and rho)
sig_inc = sigma/((1-rho^2)^0.5);

(lNormDev, ly ) = getEquiprobNormalDeviates(mu, sig_inc, numPointsY ); # thanks to Alex

#---------------------%
#Get transition matrix Q(i, j). The prob of income j in t+1
#conditional on income i in t
#---------------------%
 Q = Array(Float64,numPointsY, numPointsY);             #initialise the transition matrix


for i = 1:1:numPointsY
    for j = 1:1:numPointsY
        hiDraw = lNormDev[j+1] - (1-rho)*mu - rho * ly[i]; #highest innovation that will give us income j tomorrow
        loDraw = lNormDev[j]   - (1-rho)*mu - rho * ly[i]; #lowest  innovation that will give us income j tomorrow

        P_part1=cdf(Normal(), hiDraw/sigma );
        P_part2=cdf(Normal(), loDraw/sigma );

        Q[i,j] = P_part1 - P_part2;
    end #j

    #Each of the rows of Q should add up to 1. But
    #due to the truncation of the normal distribution they will
    #not. So we divide through by the sum of elements in the row
    Q[i, :] = Q[i, :] ./ sum(Q[i, :]);
end #i

y = exp(ly);                # Get y from log y
minInc = exp(-normBnd * sig_inc); # Get the minimum income in each year
maxInc = exp(normBnd * sig_inc);  # Get the maximum income in each year

if (y[1] < 1e-4) || (y[numPointsY] > 1e5)
    println("Combination of sigma and rho give a very high income variance. Numerical instability possible")
end


#----------------------------------------%
# Now get a matrix, T * numPointsY that holds the grid for each income in
# each year. Do likewise with minimum and maximum income
#----------------------------------------%
Ygrid = repmat(y', T, 1);
minInc = repmat([minInc]', T, 1);
maxInc = repmat([maxInc]', T, 1);

#----------------------------------------%
#Replace these arays with zeros for all years after retirement
#----------------------------------------%
if Tretire == 0         #no work (retired at birth)
   Ygrid[:, :] = 0;
   minInc[:, :] = 0;
   maxInc[:, :] = 0;
elseif (Tretire > 0) & (Tretire <=T)  #retire at some age
    Ygrid[Tretire:T, :] = 0;
    minInc[Tretire:T, :] = 0;
    maxInc[Tretire:T, :] = 0;
end

println("Income Grid Done")

return Ygrid, Q, minInc, maxInc

end
