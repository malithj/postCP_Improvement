## ---- fig.width=5, fig.align='center'------------------------------------
require(postCP)
#prepare data
sigma=1.3
#Change point estimates
bp=c(7,10)
#Obtain data from longley dataset
data = longley
#Apply postcp function
res = postcp(Employed ~ GNP + Armed.Forces,family=gaussian(),data=data,bp=c(7,10),sigma)

#Plot the results
plot.postCP(res,main="Posterior Change Point Probability Distribution")

#Apply postcp function with maxFB=TRUE to obtain marginal distribution
res = postcp(Employed ~ GNP + Armed.Forces,family=gaussian(),data=data,bp=c(7,10),sigma,maxFB=TRUE)

#Plot the results
plot.postCP(res,main="Posterior Change Point Probability Distribution")

#plot.postCP also provides support to display the response variable easily
plot.postCP(res,show.response = TRUE)

