# Paul Rodriguez translation into Julia 0.3.5 of Monica Costa Dias and Corma O'Dea "Dynamic Economics in Practice"
# using Alex code for Getting Equiprobable Shock Points from a Normal Distributions
# March 2015


# ------------------------------------------------------------------------
# DESCRIPTION
# This program solves and simulates a finite period consumption and saving
# problem. There is income that can be uncertain. The income program can be
# hardcoded by the user or can be set to follow a log normal autoregessive
# process

# ------------------------------------------------------------------------
# PREAMBLE
# Ensure that all storage spaces variables, matrices, memory, globals are
# clean from information stored in past runs

using DataArrays	  # Handle missings in the simulations
using DataFrames      # For making graphs or exporting data, is better as dataframes

using Distributions   # For drwaing distributions
using NLopt           # Basic optimization routines
using Grid            # Interpolation tools


#cd("C:\\Users\\PaulAndrés\\Dropbox\\Health and Labour Supply\\Model\\Dynamic Economics 28-29 October 2013\\final code\\v5_julia")
cd("C:\\Dropbox\\Health and Labour Supply\\Model\\Dynamic Economics 28-29 October 2013\\final code\\v5_julia")
pwd()



# Numeric tools
require(joinpath("code_numericTools","getEquiprobNormalDeviates.jl"))
require(joinpath("code_numericTools","truncate.jl"))


# Grids and other inputs
require(joinpath("code_gridsInputs","checkInputs.jl"))
require(joinpath("code_gridsInputs","getIncomeGrid.jl"))
require(joinpath("code_gridsInputs","getMinAndMaxAss.jl"))
require(joinpath("code_gridsInputs","getGrid.jl"))         #It includes tLog() and eExp()

# Programs for the solution
require(joinpath("code_solution","solveValueFunction.jl"))
require(joinpath("code_solution","objectivefunc.jl"))
require(joinpath("code_solution","utility.jl"))

# Programs for the simulation
require(joinpath("code_simulation","simNoUncer.jl"))
require(joinpath("code_simulation","simWithUncer.jl"))

# Plot stuff
require(joinpath("code_plots","plotPaths.jl"))
require(joinpath("code_plots","plots.jl"))


tic()        # start the clock

# ------------------------------------------------------------------------
# DECLARE VARIABLES AND MATRICES THAT WILL BE 'GLOBAL'
# explicitly set variables and matrices to be shared throughout the routine
# as globals

#% ------------------------------------------------------------------------
# NUMERICAL METHODS
# select solution, interpolation and integration methods

# ------------------------------------------------------------------------
# NUMERICAL CONSTANTS
# set constants needed in numerical solution and simulation

# precision parameters
#--------------------------------------%
const tol = 1e-10;                 # max allowed error
const minCons = 1e-5;              # min allowed consumption

# where to truncate the normal distributions
#--------------------------------------%
const normBnd = 3;                 #Ignore draws less than -NormalTunc*sigma and greater than normalTrunc*sigma


# information for simulations
#--------------------------------------%
const numSims = 2;                #How many individuals to simulate


# ------------------------------------------------------------------------
# THE ECONOMIC ENVIRONMENT
# Set values of structural economic parameters

const T = 50;                      	# Number of time period
const r = 0.01;                    	# Interest rate
const borrowingAllowed = 1;        # Is borrowing allowed
const isUncertainty = 1;           	# Is income uncertain?
const startA = 0;                  	# How much asset do people start life with
Tretire = 41;                		# age after which there is no income earned

beta = 0.98;                 		# Discount factor
gamma = 1.5;                 		# Coefficient of relative risk aversion
paramsUtil=[beta,gamma]

mu = 0;                     		# mean of initial log income
sigma = 0.25;                   	# variance of log income
rho = 0.75;                     	# persistency of log income
paramsInco=[mu,sigma,rho]


# Construct a vector of parameters
thetaL= [length(paramsUtil),length(paramsInco)]
theta 		= Array(Float64,sum(thetaL)) # Big vector of parameters

#Feed theta
theta[                 1:thetaL[1]]=paramsUtil; 		# 7 params
theta[       thetaL[1]+1:sum(thetaL[1:2])]=paramsInco; 	# 12 params


# ------------------------------------------------------------------------
# GRIDS
# choose dimensions, set matrices and select methods to construct grids

#The grid for assets
#--------------------------------------%
const numPointsA = 20;             # number of points in the discretised asset grid
gridMethod = "3logsteps";    # method to construct grid. One of equalsteps, logsteps, 3logsteps, 5logsteps or 10logsteps

#The grid for income shocks
#--------------------------------------%
const numPointsY = 5;           #  points in grid for income (should be 2 if hard-coded)
const uncertaintyMethod = 1;    #  =0 if we enter shocks manually, =1 if put shocks on a grid using Tauchen (1986) method
hcIncome = [0.5941 ,   0.8179  ,  1.0000  ,  1.2227  ,  1.6832];  # hard-coded shocks, used if uncertaintyMethod = 1
hcIncPDF = [0.1016 ,   0.2492  ,  0.2983  ,  0.2492  ,  0.1017];  # and respective probabilities

## Get income grid
checkInputs(paramsInco)
(Ygrid, incTransitionMrx, minInc, maxInc) = getIncomeGrid(paramsInco);


## ------------------------------------------------------------------------
# GET ASSET GRID
# populate grid for assets using 'gridMethod'
# populate matrix borrowingConstraints

(borrowCon, maxAss) = getMinAndMaxAss(borrowingAllowed, minInc, maxInc, startA);

Agrid = Array(Any,T+1); # This is a grid of log-assets shifted (so the negatives are not a problem)!! so you must be careful on recover the real A
										# every moment that you want to use it
for ixt = 1:1:T+1
	Agrid[ixt] = getGrid(borrowCon[ixt], maxAss[ixt], numPointsA, gridMethod)
end

 eExp( [Agrid[1]] ,borrowCon[1],gridMethod)

#% ------------------------------------------------------------------------
# SOLVE CONSUMER'S PROBLEM
# Get policy function and value function

(policyA1, policyC, val, EV) = solveValueFunction(theta);

	# Want to compare with Matlab version? In that way you can test how different are
	# simulations for the same policyFunc in Matlab and Julia
	# policyA1=reshape(readdlm("matlablObj\\policyA1matlab.csv", ',' ),40,20,5);

#% ------------------------------------------------------------------------
# SIMULATE CONSUMER'S PATHS
# start from initial level of assets and simulate optimal consumption and
# savings profiles over lifecycle

if isUncertainty == 0
    data1F = simNoUncer(theta,policyA1, EV, startA);  # data1F= (ypath, cpath, apath, vpath)
else
	# Get the draws, and the feed the simulation with them (sligthly different than the original code)

	# For the discreate shocks of Income case... Julia makes our life way easier!
	srand(1223424) # Seed it!
	randYpath=rand( Categorical(vec(hcIncPDF))  ,(T,numSims));

	# For the normal distributed case
	srand(1223424) # Seed it!
 	e=rand(Normal(0, sigma),(T,numSims)); # normally distributed random draws for the innovation
 	srand(234636)  # Seed it!
 	sig_inc = sigma/ ((1-rho^2)^0.5);
 	logy1   =rand(Normal(mu, sig_inc),numSims);  # a random draw for the initial income

		# Want to compare with Matlab version? In that way you can test how different are
		# simulations using the same shocks
 		e      = readdlm("matlablObj\\ematlab.csv", ',' );
 		logy1  = readdlm("matlablObj\\logy1matlab.csv", ',' );

    Ret1 = simWithUncer(theta,policyA1,EV, startA,e,logy1,randYpath);  # data1F= (ypath, cpath, apath, vpath)
end


toc();     # Stop the clock


#% ------------------------------------------------------------------------
# Export data

 # Return the simulated dataset dataset

    indi=[1:numSims]
    indix=indi
    time =fill(1, numSims)
    for i=2:T
        indix = [ indix , indi ]
        time  = [ time  , fill(i, numSims) ]
    end

    data1F=DataFrame(Float64, T*numSims,length(Ret1));
    data1F[:t]=time;
    data1F[:indiv]=indix;

    for el=1:1:length(Ret1)
        a1=Ret1[el]'
        a1r=reshape(a1[:,1:50], T*numSims)
        data1F[:,el] = a1r
    end

    nameCols=[symbol("labY");symbol("C");symbol("A");symbol("V");symbol("t")]
    nameColsF=[nameCols;symbol("Indiv")]
    names!(data1F, nameColsF)

# Averages

data1=DataFrame(Float64, T,length(data1F)-1);

for el=1:1:(length(data1F)-1)
	for t=1:1:T
		a1= data1F[data1F[:t] .== t, el]
		data1[t,el] = mean(dropna(a1))
	end
end
data1
names!(data1, nameCols)

writetable(joinpath("output","data1.csv"),data1)		# All Data
writetable(joinpath("output","dataFull1.csv"),data1F)   # Just the average

#% ------------------------------------------------------------------------
# PLOTS (but is cheaper to do the in Stata, R, etc.)
# Plots some features of the simulation and simulation
# Use a program that supports Gadfly, in these days, either JUNO, IJULIA

if 1==1  # In purspouse, Gadfly is just so slow!! But there are few good alternatives (March 2015)

	using Gadfly

	plotPaths()

	# Now plot value and policy functions
	whichyear = 20;
	plotNode1 = 3;
	plotNodeLast = numPointsA;

	plots();   #Plot Value Functions
end

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
