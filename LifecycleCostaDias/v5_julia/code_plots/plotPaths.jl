function plotPaths( )


bsample = DataFrame(age=data1[:t],Sim1=data1F[ (data1F[:Indiv].==1) ,:A],Sim2=data1F[ (data1F[:Indiv].==2) , :A])
bsample[:mtindex] = 1:size(bsample,1)     # Add an identifier
datamia=stack(bsample,[2:3])              # Reshape the data so we can use colors!

plot1=plot(datamia,x="age",y="value",color="variable",Geom.line,Guide.title("Time path of assets"))
draw(PNG(joinpath("output","images","pathAssets.png"), 24cm, 12cm), plot1)

# -------------------------------------------------------------------------

bsample = DataFrame(age=data1[:t],C1=data1F[ (data1F[:Indiv].==1) ,:C],Y1=data1F[ (data1F[:Indiv].==1) , :labY])
bsample[:mtindex] = 1:size(bsample,1)     # Add an identifier
datamia=stack(bsample,[2:3])              # Reshape the data so we can use colors!

plot2=plot(datamia,x="age",y="value",color="variable",Geom.line,Guide.title("Time path of income and consumption Individual 1"))
draw(PNG(joinpath("output","images","pathIndiv1.png"), 24cm, 12cm), plot2)

# -------------------------------------------------------------------------

bsample = DataFrame(age=data1[:t],C2=data1F[ (data1F[:Indiv].==2) ,:C],Y2=data1F[ (data1F[:Indiv].==2) , :labY])
bsample[:mtindex] = 1:size(bsample,1)     # Add an identifier
datamia=stack(bsample,[2:3])              # Reshape the data so we can use colors!

plot3=plot(datamia,x="age",y="value",color="variable",Geom.line,Guide.title("Time path of income and consumption Individual 2"))
draw(PNG(joinpath("output","images","pathIndiv2.png"), 24cm, 12cm), plot3)

end


