
function checkInputs(paramsInco)

mu    = paramsInco[1];
sigma = paramsInco[2];
rho   = paramsInco[3];
    
#Check various inputs

 println("Let's check inputs...!")


# Check that the value of rho is not greater than 1. Warn if close to 1
if (rho>0.999) || (rho <-0.999)
    if (rho>1)  || (rho < 1)     
    error("rho is greater than 1. This code solves a stationary income process. rho greater than 1 implies a non-stationary process")
    else
    warning("rho is greater than 0.99. This code solves a stationary income process. As rho gets closer to 1 the process becomes nonstationary - possibility of numerical instability")        
    end
end

#
# Check that the standard deviation is not too small
if (sigma<1e-10) 
    if (sigma<=0) 
        error("sigma is less than or equal to zero")
    else
        warning("sigma is very small and close to zero - possibility of numerical instability. Consider turning uncertainty off")
    end
end


if (borrowingAllowed != 0) && (borrowingAllowed !=1)
    error("borrowingAllowed should be either 0 or 1")
end

println("Inputs are fine :)")

end
