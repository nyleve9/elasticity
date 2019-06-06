#xx=matrix(rep(stable.stage,3),3,3)
#xx
#ss.matrix=t(xx) #flips rows to columns
#ss.matrix
#repro.matrix=matrix(rep(repro.value,3),3,3) #does the repro.matrix need to be flipped, too?
#repro.matrix

#---------making a matrix of sensitivities------------------------------------

possum<-matrix(c(0,0.616,0.616,1,0,0,0,0.7,0.7),3,3,byrow=TRUE)
possum
eigen(possum)

source("eigenall.R")
eigenall(possum)
outs=eigenall(possum)
outs

stable.stage<-(eigenall(possum)$stable.stage)
stable.stage
repro.value<-(eigenall(possum)$repro.value) 
repro.value

sens.denom<-sum(stable.stage*repro.value)
sens.denom

w<-cbind(stable.stage) #makes a vertical array
w
v<-rbind(repro.value) #makes a horizontal array
v

#to calculate sensitivities of all the matrix elements:
sensitivity<-(w%*%v)/sens.denom   #is this right?
sensitivity

#add up the relevant matrix elements to get sensitivities to vital rates
#then convert sensitivities to elasticities

#--------calculating elasticities, brute force---------------------
possum<-matrix(c(0,.616,.616,1,0,0,0,.7,.7),3,3,byrow=TRUE)
possum
eigen(possum)
possum.as<-matrix(c(0,.616,.585,1,0,0,0,.7,.665),3,3,byrow=TRUE)
possum.as
eigen(possum.as)
possum.js<-matrix(c(0,.585,.616,1,0,0,0,.665,.7),3,3,byrow=TRUE)
possum.js
eigen(possum.js)
possum.jo<-matrix(c(0,.585,.616,1,0,0,0,.7,.7),3,3,byrow=TRUE)
possum.jo
eigen(possum.jo)
possum.ao<-matrix(c(0,.616,.585,1,0,0,0,.7,.7),3,3,byrow=TRUE)
possum.ao
eigen(possum.ao)
#----------------------------------------------------
possum<-matrix(c(0,.616,.616,1,0,0,0,.7,.7),3,3,byrow=TRUE)
pmatrix<-matrix(nrow=3,ncol=3)  #makes an empty matrix
pmatrix<-possum   #stores the possum matrix in the pmatrix
pmatrix

lams=NULL
jsurv<-seq(0.05,1,by=0.01)
for(i in jsurv) {
  pmatrix[3,2]<-i
  pmatrix[1,2]<-i*0.88
  ev<-eigen(pmatrix)
  lams<-c(lams,ev$values[1])
}
lams<-as.matrix(lams)
jsurv<-as.matrix(jsurv)
plot(jsurv,lams,xlab="juvenile survival",ylab="lambda",type="l")

plams=NULL
percjs=NULL
l=NULL
js=NULL
for(i in 1:length(jsurv)) {
  l<-((lams[i+1]-lams[i])/lams[i])
  js<-((jsurv[i+1]-jsurv[i])/jsurv[i])
  plams<-c(plams,l)
  percjs<-c(percjs,js)
}
percjs<-percjs*100
plams<-plams*100
percjs
plams
plot(percjs,plams,xlab="percent change in juvenile survival",ylab="percent change in lambdas",type="l")

plams=NULL
percjs=NULL
l=NULL
js=NULL
for(i in 1:length(jsurv)) {
  l<-((lams[i]-1.209)/1.209)
  js<-((jsurv[i]-0.7)/0.7)
  plams<-c(plams,l)
  percjs<-c(percjs,js)
}
percjs<-percjs*100
plams<-plams*100
plot(percjs,plams,xlab="percent change in juvenile survival",ylab="percent change in lambdas",type="l")

jsurv<-seq(0.05,1,by=0.01)
count<-NULL
cc<-NULL
for(i in jsurv) {
  cc<-(1-i)/i
  count<-c(count,cc)
}
count<-count*100
plot(count,plams,type="l")

#----------------NAKED STINK RAT------------------------
mx1<-matrix(c(0,0.2,4,0.1),2,2)
mx2<-matrix(c(0,0.3,0.4,0.8),2,2)
mx3<-matrix(c(0.8,0.05,0.01,0.1),2,2)
# eigen(mx1)
# eigen(mx2)
# eigen(mx3)
# allmxs<-matrix(nrow=2,ncol=6)
# allmxs[1:2,1:2]<-mx1
# allmxs[1:2,3:4]<-mx2
# allmxs[1:2,5:6]<-mx3
# allmxs
# sample(allmxs,mx1,mx2,mx3)

#SimpleGrowthChooser.r
# This program allow you to simulate a 
# very simple growth process, with you providing a discrete
# list of lambda values and probabilities that each will occur

#graphics.off()# closes fig window
rm(list=ls())  # clear all variables
mx1<-matrix(c(0,0.2,4,0.1),2,2)
mx2<-matrix(c(0,0.3,0.4,0.8),2,2)
mx3<-matrix(c(0.8,0.05,0.01,0.1),2,2)

# -----INPUT PARAMETERS----------------------

maxyr = 100 # the number of years to simulate

startN = c(500,500) # starting population size vector/matrix

Reps = 5 # number of replicate runs to make  	

#---- END OF INPUTS ---------------------------
Ns = matrix(0,maxyr,Reps) # set up a matrix for the pop sizes to fill
n<-c(1,2,3)

for (jj in 1:Reps) {
  N0 = startN
  Ns[1,jj]=sum(startN) # initialize with starting pop size
  for (ii in 1:(maxyr-1)) {
    #lam_t=lams[sample(ncol(lams),size=2,replace=TRUE) # choose a random lambda value
    mx=NULL
    choose <- sample(n,1,replace=TRUE)
    if (choose == 1 ) {
          mx = mx1
    } else if (choose == 2) {
          mx = mx2
    } else {
          mx = mx3
    }           
    
    N0<-mx%*%N0
    Ns[(ii+1),jj]<-sum(N0)
    
  } # end of ii loop
} # end of jj loop

l=length(Ns)
l
lambdas<-((Ns[2:l])/(Ns[1:(l-1)]))
mean(log(lambdas))
lambdas
Ns

#perparing a matrix to use in estimating extinction risk:
Ns2=Ns 
Ns2[Ns2 >0] = 1 #set all values of Ns that are greater than one to 
alive = apply(Ns2,1,sum)
dead = 100*(1-alive/Reps)


allyrs = c(1:maxyr)
windows()

plot(allyrs,dead,type = "l",xlab="Year",ylab="Extinction CDF")

windows()
matplot(allyrs,Ns, type = "l",xlab="Year",ylab="Nt",main="Arithmetic scale")

windows()
matplot(allyrs,log(Ns), type = "l",xlab="Year",ylab="log(Nt)",main="Log scale")