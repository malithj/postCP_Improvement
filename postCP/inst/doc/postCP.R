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

## ---- fig.width=5, fig.align='center'------------------------------------
require(postCP)
data = longley
plot(data$Armed.Forces)
res = postcp(Armed.Forces ~ 1,family=gaussian(),data=data,bp=c(10),sigma=1)
plot.postCP(res,main="Posterior Change Point Probability Distribution")


## ---- fig.width=5, fig.align='center'------------------------------------

## change in the mean
data = data.frame(signal=rnorm(150) + rep(c(1, 4), each=75), position=1:150)
plot(data$signal)
res = postcp(signal ~ 1 ,family=gaussian(),data=data,bp=c(40),sigma=1)
plot.postCP(res,main="Posterior Change Point Probability Distribution")

## change in the slope
position <- 1:150
data = data.frame(signal = rep(c(1, 2), c(50, 100))*position + rnorm(150), position=position)
plot(data$signal)
res = postcp(signal ~ position ,family=gaussian(),data=data,bp=c(30),sigma=1)
plot.postCP(res,main="Posterior Change Point Probability Distribution")


## ---- fig.width=5, fig.align='center'------------------------------------

## change in the mean with number = 10000
number = 10000
data = data.frame(signal=rnorm(number) + rep(c(2, 5), each=number/2), position=1:number)
plot(data$signal)
res = postcp(signal ~ 1 ,family=gaussian(),data=data,bp=c(3000),sigma=1)
plot.postCP(res,main="Posterior Change Point Probability Distribution")


