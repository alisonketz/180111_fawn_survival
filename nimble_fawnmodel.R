#Dan Walsh
#Date: 2/8/18

rm(list=ls())

# setwd('F:/171206_Survival')

library(rjags)
library(dclone)
library(doParallel)
library(foreach)
library(coda)
library(Hmisc)
library(lubridate)
library(xlsx)
library(gridExtra)
library(xtable)

###
### Preliminaries
###

setwd('/home/aketz/Documents/Survival/180111_fawn_survival')
load('jags.data.Rdata')

####Nimble Test
library(nimble)
library(ggplot2)
library(beepr)
library(Matrix)
library(coda)

simcode <- nimbleCode({
  
  # Priors
  beta0 ~ dnorm(0, 1.0E-6) # intercept log hazard rate
  
  # Likelihood for the total hazard
  for (j in 1:records) {
    for (k in left[j]:(right[j]-1)) {
      UCH[j,k] <- exp(beta0 + timeuse[k])
    }
    
    # total prob of surviving
    SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
    censor[j] ~ dbern(SLR[j])      
  }
  
  timeuse[1:nT] ~ dcar_normal(adj[1:nN],weights=weights[1:nN],num[1:nT],tau, zero_mean=1)
  
  tau ~ dgamma(1,1)
  
  #Derived parameters
  for (t in 1:nT){
    llambda.out[t]<-beta0 + timeuse[t] #CWD-
    UCH0[t]<-exp(llambda.out[t])
    CH0[t]<-sum(UCH0[1:t])
    S0[t]<-exp(-CH0[t])
    
  }
})

#create num vector
num<-c(1,rep(2,jags.data$nT-2),1)

#create adj vector
temp<-as.matrix(bandSparse(n=jags.data$nT,k=c(1),symmetric = T))
temp2<-matrix(0,jags.data$nT,jags.data$nT)
for(i in 1:nrow(temp2)){
  temp2[i,] <-  temp[i,]*c(1:jags.data$nT)
}
adj<-t(temp2)[which(t(temp2)!=0)]
adj
simData <- list(censor=jags.data$censor) 
simConsts <- list(x5=jags.data$x5, nT=jags.data$nT, records=jags.data$records, 
                  right=jags.data$right, left=jags.data$left,adj=adj, num=num, nN=length(adj),weights=rep(1,length(adj)))
simInits <- list("beta0"=rnorm(1,0,1),timeuse=rep(1,jags.data$nT))

simtest<-nimbleModel(code= simcode, name="fawn", constants = simConsts,data = simData, inits = simInits)

###############################################################################RW config

# reps<-200000
# bin<-0.1*reps
# nchains<-3
# 
# ##One-line MCMC call
# aaa<-proc.time()
# mcmcout<-nimbleMCMC(model=simtest, nchains=nchains, nburnin = bin,
#                     niter=reps, summary=TRUE,
#                     monitors=c("beta0","tau","S0","CH0")) 
# proc.time()-aaa
# 
# beep(sound=4)
# 
# 
# out<-mcmc.list(lapply(mcmcout$samples, mcmc))
# traceplot(out[,"beta0"],ylab="beta0")
# traceplot(out[,"tau"],ylab="tau")
# densplot(out[,"tau"])

###############################################################################Manual config

reps<-200000
bin<-0.1*reps
nchains<-3

mcmcconf<-configureMCMC(simtest,
                        monitors=c("beta0","tau","S0","CH0"),print=T)
mcmcconf$setSamplers(2:3)
mcmcconf$addSampler(target="beta0",type="slice")
mcmcconf$printSamplers()

Rmcmc<-buildMCMC(mcmcconf)
Cmcmc<-compileNimble(Rmcmc)

aa<-proc.time()
mcmcout2<-runMCMC(Cmcmc, nchains=nchains, nburnin = bin,
                  niter=reps, summary=TRUE)
proc.time()-aa

for(i in 1:10){
  beep(sound=4)
}

out<-mcmc.list(lapply(mcmcout2$samples, mcmc))
traceplot(out[,"beta0"],ylab="beta0")
traceplot(out[,"tau"],ylab="tau")

fit.sum=summary(out)[[1]]
fit.quant=summary(out)[[2]]
tail(fit.sum)
tail(fit.quant)

save.image("nimble.fawn.Rdata")


gelman.diag(out)



######################################################################################################################3
###
### Output
###
######################################################################################################################3

###
### Notable Days
###
#day cut-off for harvest 

start.day=yday("2017-05-16")
non.harvest.survival.end=yday('2017-09-16')-start.day
bow.start.day=yday('2017-09-16')-start.day+1
bow.end.day = yday('2018-01-07')+365-start.day+1
gun.start.day=yday('2017-11-18')-start.day+1
gun.end.day = yday('2017-11-26')-start.day+1
gun.holiday.start.day = yday('2017-12-24')-start.day+1
gun.holiday.end.day = yday('2018-01-01')+365-start.day+1
study.end=yday('2018-01-07')+365-start.day+1
n.days=study.end





# 
# 
# n.sum=dim(fit.sum)[1]
# round(fit.sum[(n.sum-188):n.sum,1])
# which(N.data.fit[,4]==1)
# sum(N.data.fit[,4])
# z.out=as.numeric(fit.sum[(n.sum-188):n.sum,1]>.2)
# sum(z.out)
# 
# which(z.out==1)


# The Colorblind palette with grey:
cbPalette <- c( "#0072B2", "#D55E00", "#CC79A7","#999999", "#E69F00", "#56B4E9")

###
### plotting single survival curve
###

Days=rep(1:n.days)
Hazard=fit.sum[1:n.days,1]
CH0.l=fit.quant[1:n.days,1]
CH0.u=fit.quant[1:n.days,5]

out.ch0= data.frame(Days,Hazard,CH0.l,CH0.u)

ggplot(data = out.ch0, aes(x = Days,y=Hazard))+
    geom_ribbon(aes(ymin=CH0.l,ymax=CH0.u),alpha=.1,show.legend=NA,linetype=0)+
    geom_line()+theme_bw()


Survival=fit.sum[(n.days+1):(2*n.days),1]
Lower=fit.quant[(n.days+1):(2*n.days),1]
Upper=fit.quant[(n.days+1):(2*n.days),5]

out.df=data.frame(cbind(Days,Survival,Lower,Upper))


ggplot(data = out.df, aes(x = Days,y=Survival))+
    geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=.1,show.legend=NA,linetype=0)+
    geom_line()+theme_bw()

pdf("Survival_fawn_v3.pdf", width=8, height =6)
ggplot(data =out.df,aes(x = Days,y=Survival))+
    annotate('rect',xmax=gun.holiday.end.day, xmin=gun.holiday.start.day, ymax=Inf, ymin=-Inf, linetype=2,alpha=.5,fill="grey90")+
    annotate('rect',xmax=gun.end.day, xmin=gun.start.day, ymax=Inf, ymin=-Inf, linetype=2,alpha=.5,fill="grey90")+
    geom_line()+theme_bw()+
    geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=.1,show.legend=NA,linetype=0)+
    ggtitle("Fawn Survival")+xlab("Time (Days)")+ylab("Survival Probability")+
    geom_vline(aes(xintercept=bow.start.day),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.start.day),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.end.day),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.holiday.start.day),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.holiday.end.day),linetype=2,color="grey50")+
    scale_x_discrete(limit = c(1,51,101,151,gun.start.day,gun.holiday.start.day),labels = c("May 16","July 5","Aug 24","Oct 13","Nov 18","Dec 24","Jan 7"))
dev.off()

###
### Hazard Plot
###

pdf("Hazard_fawn_v3.pdf", width=8, height =6)
ggplot(data =out.ch0,aes(x = Days,y=Hazard))+
    annotate('rect',xmax=gun.holiday.end.day, xmin=gun.holiday.start.day, ymax=Inf, ymin=-Inf, linetype=2,alpha=.5,fill="grey90")+
    annotate('rect',xmax=gun.end.day, xmin=gun.start.day, ymax=Inf, ymin=-Inf, linetype=2,alpha=.5,fill="grey90")+
    geom_line()+theme_bw()+
    geom_ribbon(aes(ymin=CH0.l,ymax=CH0.u),alpha=.1,show.legend=NA,linetype=0)+
    ggtitle("Fawn Hazard")+xlab("Time (Days)")+ylab("Cumulative Hazard")+
    geom_vline(aes(xintercept=bow.start.day),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.start.day),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.end.day),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.holiday.start.day),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.holiday.end.day),linetype=2,color="grey50")+
    scale_x_discrete(limit = c(1,51,101,151,gun.start.day,gun.holiday.start.day),labels = c("May 16","July 5","Aug 24","Oct 13","Nov 18","Dec 24","Jan 7"))
dev.off()



###
### Hazard analysis coefficients
### 

output=cbind(fit.sum[(2*n.days+1):(2*n.days+2),1:2],fit.quant[(2*n.days+1):(2*n.days+2),c(1,5)])
row.names(output)=c("Intercept",expression(tau))
output=round(output,2)
output

write.csv(output,file="Hazard_table_fawn.csv",row.names = TRUE)

output2=cbind(c("Intercept","$//tau$"),output)

print.xtable(xtable(output2,caption = 'Coefficients of the hazard rate analysis where the baseline hazard shows the daily hazard intercept, and additive term for Gun Harvests. The first column denotes the mean of the posterior distribution with the standard deviation (SD), and equal-tailed Bayesian credible intervals in the right two columns.',
                    align = c('c','c','c','c','c','c'),label="tab:hazard_fawn",digits=c(0,0,2,2,2,2)),
             sanitize.text.function = function(x) {x},
             hline.after=c(0,2),
             include.rownames=FALSE,
             rotate.colnames = FALSE,
             caption.placement = "top")


###################################################################################################
###
### the "annual" survival estimates (i.e., the mean/median and 95%CI for surviving to the end of the analysis period) broken down by harvest types
###
###################################################################################################

###
### Table with survival before/after bow and gun hunts
###
pre.end=non.harvest.survival.end
bow.end=gun.start.day-1
gun1.end=gun.end.day
gun2.end=gun.holiday.end.day

fit.sum[c(study.end+bow.end),1]
pre.bow.harvest=c(fit.sum[study.end+pre.end,1:2],fit.quant[study.end+pre.end,c(1,5)])
pre.bow.harvest=round(pre.bow.harvest,3)
pre.bow.harvest
pre.bow.harvest=data.frame(rbind(pre.bow.harvest))
names(pre.bow.harvest)=c("Mean","SD","0.025","0.975")
names(pre.bow.harvest)

post.gun1.harvest=cbind(fit.sum[study.end+gun1.end,1],fit.sum[study.end+gun1.end,2],fit.quant[study.end+gun1.end,1],fit.quant[study.end+gun1.end,5])
post.gun1.harvest=round(post.gun1.harvest,3)
post.gun1.harvest
post.gun1.harvest=data.frame(post.gun1.harvest)
names(post.gun1.harvest)=c("Mean","SD","0.025","0.975")


post.gun2.harvest=cbind(fit.sum[study.end+gun2.end,1],fit.sum[study.end+gun2.end,2],fit.quant[study.end+gun2.end,1],fit.quant[study.end+gun2.end,5])
post.gun2.harvest=round(post.gun2.harvest,3)
post.gun2.harvest
post.gun2.harvest=data.frame(post.gun2.harvest)
names(post.gun2.harvest)=c("Mean","SD","0.025","0.975")



Period=c("Pre-Gun Harvest","Gun","Holiday Gun")

survival_all_sum=data.frame(cbind(Period,rbind(pre.bow.harvest,post.gun1.harvest,post.gun2.harvest)))

survival_all_sum
names(survival_all_sum)=c("Period","Mean","SD","0.025","0.975")

write.csv(survival_all_sum,"Survival_fawn.csv",row.names=FALSE)

print.xtable(xtable(survival_all_sum,caption = 'The mean of the posterior distribution (black line) of the daily cumulative probability of survival of fawns after birth during 2017 of white-tailed deer in southwest Wisconsin, along with 95\\% equal-tailed Bayesian credible intervals (gray shaded region).', 
                    align = c('c','c','c','c','c','c'),label="tab:survival_fawn",digits=c(0,0,2,2,2,2)),
             sanitize.text.function = function(x) {x},
             hline.after=c(-1,-1,0,3,3),
             include.rownames=FALSE,
             rotate.colnames = FALSE,
             caption.placement = "top")  


#saveworking directory
save.image("fawn_survival_v3.Rdata")

