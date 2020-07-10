# The goal of this analysis is to explore the following
# question. How much of the expected survival probability
# needs to be explained by the biomarker for it to be useful
# for adpative dosing. Adaptive dosing here relates to dosing
# being turned on/off based on the value of a biomarker. 

# The model we will use is an L-V model with a drug component.

# Model
require(deSolve)
LV <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    dD <- 0                                # Drug
    dS <- gS*S*(1-(S+Krs*R)/K1) - d*S*D    # Drug Sensitive
    dR <- gR*R*(1-(R+Ksr*S)/K2)            # Drug Resistant
    dES <- -(a1*D + a2*(S+R))*ES           # Expected Survival Probability  
    res <- c(dD, dS, dR, dES)
    list(res)
  })
}

# 1st simulation with no drug - we want S to be the dominant 
# population and we allow the survival to depend on the 
# biomarker

parms<-c(gS = 0.0033, gR = 0.0017, 
         Krs = 1, Ksr = 1, K1 = 1, K2 = 1, 
         a1 = 0, a2 = 0.005,
         d = 0.0133)
# a2 is chosen to give a median OS of 0.5 with no treatment


# initial conditions
ini<-c(D = 0, S=0.75*0.73, R = 0.75*0.27, ES = 1)

# times
t<- seq(0,7*365,by=1) 

# simulate
out <- data.frame(ode(y = ini,times = t,func = LV, 
                      parms = parms))

# species
plot(out$time,out$S+out$R,type="l",
     xlab="Time (Days)",ylim=c(0,2),ylab="species")
lines(out$time,out$S,col="green")
lines(out$time,out$R,col="red")
control<-out

# survival probability
plot(out$time/365,out$ES,type="l",
     xlab="Time (Years)", ylab="Survival Probability",
     ylim=c(0,1))
out$time[out$ES<=0.5][1]

# now let's add drug - continuous dosing
ini<-c(D = 1, S=0.75*0.73, R = 0.75*0.27, ES = 1)
out <- data.frame(ode(y = ini,times = t,func = LV, 
                      parms = parms))
continuous<-out

# now adaptive

# so we have a treatment effect which is captured
# by the biomarker, what would the value of a1
# need to be to give a similar survival curve

# let's use 30% reduction
threshold1<-0.75-0.75*0.40
threshold2<-0.75

ini<-c(D = 1, S=0.75*0.73, R = 0.75*0.27, ES = 1)
out <- data.frame(ode(y = ini,times = t,func = LV, 
                      parms = parms))
# now turn off drug when threshold reached
idx<-which((out$S+out$R)<threshold1)[1]
s1<-out[out$time<out$time[idx],]
ini<-c(D = 0, S=out$S[idx], R = out$R[idx], ES = out$ES[idx])
out <- data.frame(ode(y = ini,times = t,func = LV, 
                      parms = parms))
idx<-which((out$S+out$R)>threshold2)[1]
flag<-0
# keep cycling until we no longer see a 30% reduction at which
# point patient stays on treatment
while (flag==0){
  if (is.na(idx)==F){
    dummy<-out[out$time<out$time[idx],]
    dummy$time<-dummy$time+max(s1$time)
    s1<-rbind(s1,dummy)
    ini<-c(D = 1, S=out$S[idx], R = out$R[idx], ES = out$ES[idx])
    out <- data.frame(ode(y = ini,times = t,func = LV, 
                          parms = parms))
    idx<-which((out$S+out$R)<threshold1)[1]
  }else{
    out$time<-out$time+max(s1$time)
    s1<-rbind(s1,out)
    flag<-1
  }
  if(is.na(idx)==F){
    dummy<-out[out$time<out$time[idx],]
    dummy$time<-dummy$time+max(s1$time)
    s1<-rbind(s1,dummy)
    ini<-c(D = 0, S=out$S[idx], R = out$R[idx], ES = out$ES[idx])
    out <- data.frame(ode(y = ini,times = t,func = LV, 
                          parms = parms))
    idx<-which((out$S+out$R)>threshold2)[1]
    
  }else{
    out$time<-out$time+max(s1$time)
    s1<-rbind(s1,out)
    flag<-1
  }
}

plot(control$time/365,control$S+control$R,
     xlab="Time (Years)",ylab="Total Tumour Burden",
     ylim=c(0,1),type="l",col="red")
lines(continuous$time/365,continuous$S+continuous$R,col="blue")
lines(s1$time/365,s1$S+s1$R,col="green")
mtext("No therapy",col="red",at=5,padj=2.2)
mtext("Continuous",col="blue",at=4,padj=5)
mtext("Adaptive",col="green",at=6,padj=10)

# We can see that time needed to reach the value 1 is clearly
# longer with adaptive versus continuous but what happens if we
# account for the fact that benefit is conditional on you 
# surviving up until then...

# plot survival
plot(control$time/365,control$ES, type="l",col="red",
     xlab="Time (Years)", ylab="Survival Probability")
lines(continuous$time/365,continuous$ES,col="blue")
# new drug gives survival benefit
lines(s1$time/365,s1$ES,col="green")
# Oh dear! as you'd expect its not a good idea when you take
# into account time needed to survive to benefit!
mtext("No therapy",col="red",at=5,padj=2)
mtext("Continuous",col="blue",at=5,padj=3.5)
mtext("Adaptive",col="green",at=5,padj=5)

# playing around with the threshold is interesting
# shows that adaptive therapy under this parameter
# set shows no advantage

# struggled to find any parameter set!