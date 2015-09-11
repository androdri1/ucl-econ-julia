function getGrid(minongrid, maxongrid, GridPoints, method)

#NB: This code is a modification from Chris Carroll of John's Hopkins University
# version. Be careful as Chris' code generate a nice grid that you can actually
# see. This code only generates the relevant domain, the actual action is done
# by the two functions below. This is due to the current limitations of interpolation
# packages in Julia.
#
#His website is here http://www.econ2.jhu.edu/people/ccarroll/ and his
#notes on solving dynamic models are very good - they can be found
#here: http://www.econ2.jhu.edu/people/ccarroll/SolvingMicroDSOPs.pdf


#% ------------------------------------------------------------------------
# Explanation

# We need to a grid from a to b. A basic approach could involve spacing out
# the point equally. The following line of code would achieve this:
#     grid= linspace(a, b, GridPoints);

# If we want to space them out so that the growth rate of the distance
# between spaces is equal we achieve this by spacing out the logs of grid
# equally. We could achieve this by the following line of code:
#     grid= exp(linspace(log(a), log(b), GridPoints));

# This is problematic if a<=0. The approach taken here is to get a grid
# from log(1) to log(b - a + 1), exponentiate and then subrtact 1
# so that we have a grid from 0 to b-a. If we have add a to each point we
# get the desired result - a log-spaced grid from a to b



span = maxongrid - minongrid;     # b - a



if method=="equalsteps"
  ste=(span)/(numPointsA-1)
  grid= minongrid:1:maxongrid;
elseif method== "logsteps"
  ste=(log(1+span))/(numPointsA-1)
  maxi=ste*(numPointsA-1)
  grid = log(1):ste:maxi
elseif method=="2logsteps"
  ste=(  log(1+log(1+span)) )/(numPointsA-1)
  maxi=ste*(numPointsA-1)
  grid = log(1):ste:maxi
elseif method=="3logsteps"
  ste=(  log(1+log(1+log(1+span))) )/(numPointsA-1)
  maxi=ste*(numPointsA-1)
  grid = log(1):ste:maxi
else
    error("Error in getgrid. You have entered an invalid method for choosing the distance between grid points. Method must be one of equalsteps, logsteps, 3logsteps, 5logsteps or 10logsteps.");
end

return grid;

end

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# In order to know which is the appropiate value of "A" to use in the grid. Use it for
# introducing a value in an interpolation, for instance
function tLog(a,lb,gridMethod)
  if gridMethod=="equalsteps"
      return a
  elseif gridMethod== "logsteps"
    return log(a+1-lb)
  elseif gridMethod=="2logsteps"
    return log(1+log(a+1-lb));
  elseif gridMethod=="3logsteps"
    return log(1+log(1+log(a+1-lb)));
  end
end

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# In order to know which is value of "A" that is implied by a position in the grid.

function eExp(loga,lb,gridMethod)
  if gridMethod=="equalsteps"
      return loga
  elseif gridMethod== "logsteps"
    return exp( loga )-1+lb
  elseif gridMethod=="2logsteps"
    return exp(exp( loga )-1)-1+lb;
  elseif gridMethod=="3logsteps"
    return exp(exp(exp( loga )-1)-1)-1+lb;
  end

end

