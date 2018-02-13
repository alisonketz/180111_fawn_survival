###
### Survival Analysis 
### A.Ketz 12/6/2017
###
###
###

###
### Preliminaries
###

rm(list=ls())

# setwd('F:/171206_Survival')

setwd('/home/aketz/Documents/Survival/180111_fawn_survival')


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
### Load Previous Run
### 

load('fawn_survival_v2.Rdata')

###
### Load Data and fix column names
###

#load Access Database into R using Hmisc package with the function mdb.get (must use this on Linux)
# database =mdb.get('~/Documents/Data/SWDPPdeerDB.MDB')

d.cap=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Adult_Capture_2017_2018")
d.fawncap=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Fawn Capture")
d.mort=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Mortalities")
d.cens=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Censor")
d.cwd=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "RAMALT")
d.post.cwd=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Postmortem_CWD")


names(d.cap)=tolower(gsub("[[:punct:]]","",names(d.cap)))
names(d.fawncap)=tolower(gsub("[[:punct:]]","",names(d.fawncap)))
names(d.mort)=c("InDB","Lowtag","Technicians","DataRecorder","Lefteartag","Lefteartag","Righteartag","Righteartag","Collar","Collar","Sex","Captured","MORTAlertdate","CollarFound", "EstMortDate","Comment1","PropertyType","Landtype","County","Property", "Photos","Incisor","Lymphnodes","CWDSuspect","FieldNecropsy","LabNecropsy","SecondaryFateSheet","Cause1", "Cause1wt","Cause2","Cause2wt","Cause3","Cause3wt","Comments2","HunterName","Huntingseason","Weapon","Lat","Long","CustomerID","ConfidenceinCause","NumberAntlerPoints","Township","ApproximateLocationofHarvest","NumberDeerInGroup","Howfarawaywasdeerwhenshot", "Emailaddress","CWDsampleID","Legmeasurment","Girth","Dressedweight","Hunterphone","HunterAddress","Frequency")
names(d.cens)=c("Lowtag", "Lefteartag","Righteartag","Sex","AgeAtCapture","Latitude","Longitude","County","Censorcause","CensWt","Othercause1","Othercause1wt","Othercause2","Othercause2wt","Cencomments","CaptureDate","CollarRecovered","PropertyCollarRecovered","Landtype","Frequency","CensorDate","Technicians","CollarNumber")
names(d.cwd)=tolower(gsub('[[:punct:]]',"",names(d.cwd)))
names(d.post.cwd)=tolower(gsub('[[:punct:]]',"",names(d.post.cwd)))

for(i in 1:length(names(d.mort))){names(d.mort)[i]=paste(tolower(substr(names(d.mort)[i], 1, 1)),substr(names(d.mort)[i], 2, nchar(names(d.mort)[i])),sep="")}
for(i in 1:length(names(d.cens))){names(d.cens)[i]=paste(tolower(substr(names(d.cens)[i], 1, 1)),substr(names(d.cens)[i], 2, nchar(names(d.cens)[i])),sep="")}

### Double check loading correctly

head(d.mort)
head(d.cens)
head(d.fawncap)
head(d.cwd)

###
### Number individuals within each dataframe
###

n.cap = dim(d.cap)[1]
n.cens = dim(d.cens)[1]
n.fawncap=dim(d.fawncap)[1]
n.mort = dim(d.mort)[1]
n.cwdtest=dim(d.cwd)[1]
n.postcwd=dim(d.post.cwd)[1]

#change class of CWD tests to integer for d.cwd
class(d.cwd$lowtag)="integer"

#add CWD status to capture dataframe
d.fawncap$cwdstatus=rep(NA,n.fawncap)
for(i in 1:n.fawncap){
    for(j in 1:n.cwdtest){
        if(d.fawncap$lowtag[i]==d.cwd$lowtag[j]){
            d.fawncap$cwdstatus[i]=as.character(d.cwd$result[j])
        }
    }
    if(d.fawncap$cwdstatus[i]=='ISF' | is.na(d.fawncap$cwdstatus[i]))d.fawncap$cwdstatus[i]='Negative'
}
d.fawncap$cwdstatus=as.factor(d.fawncap$cwdstatus)
d.fawncap$cwdstatus=as.numeric(d.fawncap$cwdstatus)-1
d.fawncap$cwdstatus
### 
### Dates
### Initial date of captures, ie. start of the study
### 01/09/2017
### formatting date vectors
###


d.cap$date = as.character(d.cap$date)
d.cap$date=strtrim(d.cap$date,nchar(d.cap$date)-9)
d.cap$date = as.Date(d.cap$date,format=c('%m/%d/%y'))
d.cap$dateentered = as.character(d.cap$dateentered)
d.cap$dateentered=strtrim(d.cap$dateentered,nchar(d.cap$dateentered)-9)
d.cap$dateentered = as.Date(d.cap$dateentered,format=c('%m/%d/%y'))

d.fawncap$date = as.character(d.fawncap$date)
d.fawncap$date=strtrim(d.fawncap$date,nchar(d.fawncap$date)-9)
d.fawncap$date = as.Date(d.fawncap$date,format=c('%m/%d/%y'))
d.fawncap$dateentered = as.character(d.fawncap$dateentered)
d.fawncap$dateentered=strtrim(d.fawncap$dateentered,nchar(d.fawncap$dateentered)-9)
d.fawncap$dateentered = as.Date(d.fawncap$dateentered,format=c('%m/%d/%y'))

d.mort$captured = as.character(d.mort$captured)
d.mort$captured=strtrim(d.mort$captured,nchar(d.mort$captured)-9)
d.mort$captured = as.Date(d.mort$captured,format=c('%m/%d/%y'))
d.mort$mORTAlertdate = as.character(d.mort$mORTAlertdate)
d.mort$mORTAlertdate[d.mort$mORTAlertdate==""]=NA
d.mort$mORTAlertdate[!is.na(d.mort$mORTAlertdate)]=strtrim(d.mort$mORTAlertdate[!is.na(d.mort$mORTAlertdate)],nchar(d.mort$mORTAlertdate[!is.na(d.mort$mORTAlertdate)])-9)
d.mort$mORTAlertdate = as.Date(d.mort$mORTAlertdate,format=c('%m/%d/%y'))
d.mort$estMortDate = as.character(d.mort$estMortDate)
d.mort$estMortDate[d.mort$estMortDate==""]=NA
d.mort$estMortDate[!is.na(d.mort$estMortDate)]=strtrim(d.mort$estMortDate[!is.na(d.mort$estMortDate)],nchar(d.mort$estMortDate[!is.na(d.mort$estMortDate)])-9)
d.mort$estMortDate = as.Date(d.mort$estMortDate,format=c('%m/%d/%y'))

d.mort$collarFound= as.character(d.mort$collarFound)
d.mort$collarFound[d.mort$collarFound!=""]=strtrim(d.mort$collarFound[d.mort$collarFound!=""],nchar(d.mort$collarFound[d.mort$collarFound!=""])-9)
d.mort$collarFound= as.Date(d.mort$collarFound,format=c('%m/%d/%y'))

d.cens$censorDate = as.character(d.cens$censorDate)
d.cens$censorDate=strtrim(d.cens$censorDate,nchar(d.cens$censorDate)-9)
d.cens$censorDate = as.Date(d.cens$censorDate,format=c('%m/%d/%y'))

d.cens$captureDate = as.character(d.cens$captureDate)
d.cens$captureDate=strtrim(d.cens$captureDate,nchar(d.cens$captureDate)-9)
d.cens$captureDate = as.Date(d.cens$captureDate,format=c('%m/%d/%y'))

###
### removing Adults from the mortality data frame
###

ad.mort.rm=c()
for(i in 1:n.cap){
    ad.mort.rm=c(ad.mort.rm,which(d.mort$lowtag==d.cap$lowtag[i]))
}
unique(setdiff(d.mort$lowtag,d.cap$lowtag))
d.mort[-ad.mort.rm,]


d.fawnmort=d.mort[-ad.mort.rm,]

n.fawnmort=dim(d.fawnmort)[1]


setdiff(d.fawnmort$lowtag,d.fawncap$lowtag)

which(d.fawncap$lowtag==5828)
which(d.cap$lowtag==5828)
which(d.cens$lowtag==5828)
rm.mort=which(d.mort$lowtag==5828)
d.mort=d.mort[-rm.mort,]

class(d.fawncap$lowtag)

###
### initial format of data - rows = individuals, cols = (e_i,r_i,s_i,censor)
###


#start of the study
start.day=yday("2017-05-16")


#making dates discrete for fitting model

#e
d.fawncap$dayDis = (yday(d.fawncap$date) - start.day + 1)

#r
d.fawnmort$mORTAlertdateDis=(yday(d.fawnmort$mORTAlertdate) - start.day + 1)
d.fawnmort$estMortDateDis = (yday(d.fawnmort$estMortDate) - start.day + 1)

#s
d.fawnmort$collarFoundDis=(yday(d.fawnmort$collarFound) - start.day + 1)


yday('2017-12-31')

#day cut-off for harvest 
non.harvest.survival.end=yday('2017-09-16')-start.day
bow.start.day=yday('2017-09-16')-start.day+1
bow.end.day = yday('2018-01-07')+365-start.day+1
gun.start.day=yday('2017-11-18')-start.day+1
gun.end.day = yday('2017-11-26')-start.day+1
gun.holiday.start.day = yday('2017-12-24')-start.day+1
gun.holiday.end.day = yday('2018-01-01')+365-start.day+1
study.end=yday('2018-01-07')+365-start.day+1

d.fawnmort$bow=0
d.fawnmort$bow[which(d.fawnmort$weapon!="Rifle" & d.fawnmort$cause1 == "Hunter harvest")[-1]]=1
         
# likely gun harvest
which(d.fawnmort$weapon!="Rifle" & d.fawnmort$cause1 == "Hunter harvest")

#setting harvest.season for derived parameters
n.days=study.end

# non.harvest.survival.end
# bow.start.day
# bow.end.week
# gun.start.day
# gun.end.week
# gun.holiday.start.day
# gun.holiday.end.week


bow.harvest.hazard=bow.harvest.haz = c(rep(0,non.harvest.survival.end),rep(1,bow.end.day-bow.start.day+1))
bow.harvest.hazard
length(bow.harvest.haz)

gun.harvest.haz=c(rep(0,gun.start.day-1),rep(1,gun.end.day-gun.start.day+1),rep(0,gun.holiday.start.day-gun.end.day),rep(1,gun.holiday.end.day-gun.holiday.start.day),rep(0,study.end-gun.holiday.end.day))
length(gun.harvest.haz)
gun.harvest.haz

###
### fawn # 6290 is censored but has no sex associated with it!
### just assigning it to female, although should double check, just to get model to run
###

d.fawncap$sex[52]='Female'
d.fawncap$sex=as.factor(as.character(d.fawncap$sex))


#########################################################################################################################3
###
### Formatting/combining of data
###
###
#########################################################################################################################3


###
### N.temp[,1] = e
### N.temp[,2] =r 
### N.temp[,3] = s
### N.temp[,4] = censor
### N.temp[,5] = CWD status
### N.temp[,6] = body condition
### N.temp[,7] = gun harvest
### N.temp[,8] = sex
### N.temp[,9] = bow harvest
###

#initialize matrix of captures/recaptures/morts
N.temp = matrix(NA,nr=n.fawncap,nc = 10)

#initialize with all animals alive
N.temp[,4] = 1

#initialize with all CWD neg
N.temp[,5] = 0

#initialize with gun harvest mortality=0
N.temp[,7] = 0

#initialize with sex=0 males sex=1 is females
N.temp[,8] = 0

#initialize for bow hunt
N.temp[,9] = 0


#Filling N.temp with all times
for(i in 1:n.fawncap){
    N.temp[i,1] = d.fawncap$dayDis[i]#e
    N.temp[i,5] = d.fawncap$cwdstatus[i]
    N.temp[i,6] = NA
    N.temp[i,8] = as.numeric(d.fawncap$sex[i])-1
    N.temp[i,10] = d.fawncap$lowtag[i]
    for(j in 1:n.fawnmort){
        if(d.fawncap$lowtag[i]==d.fawnmort$lowtag[j]){
            N.temp[i,2] = d.fawnmort$mORTAlertdateDis[j]#r
            N.temp[i,3] = d.fawnmort$estMortDateDis[j]#s
            N.temp[i,4] = 0
            if(d.fawnmort$cause1[j]=="Hunter harvest")N.temp[i,7]=1
            N.temp[i,9] = d.fawnmort$bow[j]
        }
        if(N.temp[i,4]==0){
            if(is.na(N.temp[i,2])){
                N.temp[i,2]=d.fawnmort$estMortDateDis[j]
            }
            if(is.na(N.temp[i,3])){N.temp[i,3]=study.end}#set to the max week, because that's the censor week of last known alive
            if(N.temp[i,2]>N.temp[i,3]){N.temp[i,2]=N.temp[i,3]}
        }
    }
}
N.temp

n.temp=dim(N.temp)[1]

###
### censoring individuals that were killed at capture
###
check.censor=c()
for(i in 1:n.fawnmort){
    check.censor=c(check.censor,which(d.cens$lowtag==d.fawnmort$lowtag[i]))
}

check.censor

#no censor morts caused by capture
#censor morts caused by capture
# mort.killcap=which(d.fawnmort$cause1=="Capture Related")
# mort.killcap
# d.fawnmort[mort.killcap,]
# cens.id=d.fawnmort$lowtag[mort.killcap]
# 
# for(i in 1:n.temp){
#     for(j in cens.id){
#         if(N.temp[i,10]==j){
#             N.temp[i,2:4]=c(study.end,study.end,1)
#         }
#     }
# }



###
### Fill in right censor
###

for(i in 1:(n.fawncap)){
    if(is.na(N.temp[i,2]))N.temp[i,2]=study.end
    if(is.na(N.temp[i,3]))N.temp[i,3]=study.end
}
tail(N.temp)


for(i in 1:n.cens){
    for(j in 1:n.temp){
        if(d.cens$lowtag[i]==N.temp[j,10]){
            cat(N.temp[j,],'\n')
        }
    }
}


#subtracting 1 day from last alive column
N.temp[N.temp[,2]==N.temp[,3] & N.temp[,2]!=study.end,2]=N.temp[N.temp[,2]==N.temp[,3] & N.temp[,2]!=study.end,2]-1



###
### Post-Mortem CWD test
###

d.post.cwd$cwdresult[d.post.cwd$cwdresult=='Pending']=NA
d.post.cwd$cwdresult=as.character(d.post.cwd$cwdresult)
d.post.cwd$cwdresult[d.post.cwd$cwdresult=='']=NA
d.post.cwd$cwdresult[d.post.cwd$cwdresult=='Negative']
d.post.cwd$cwdresult[d.post.cwd$cwdresult=='Postive']=1
d.post.cwd$cwdresult=as.numeric(as.factor(d.post.cwd$cwdresult))-1

###
### removing adults from the post mort tests
###

fawn.cwd.keep=c()
for(l in 1:n.fawncap){
    fawn.cwd.keep=c(fawn.cwd.keep,which(d.post.cwd$lowtag==d.fawncap$lowtag[l]))
}
d.post.cwd=d.post.cwd[fawn.cwd.keep,]

n.postcwd=dim(d.post.cwd)[1]


N.temp=cbind(N.temp,rep(NA,n.temp))

for(i in 1:n.postcwd){
    for(j in 1:n.temp){
         if(d.post.cwd$lowtag[i]==N.temp[j,10]){
             N.temp[j,11]=d.post.cwd$cwdresult[i]
         }
    }
}
N.temp[which(N.temp[,5]==1),11]=1


###
### Format data matrix to fit into jags
###

N.data.fit = matrix(NA,nr=n.temp+n.fawnmort-1,ncol=9)
-16
indx = 1
for(i in 1:n.temp){
    if(N.temp[i,4]==1){
        N.data.fit[indx,] = c(N.temp[i,1:2],1,N.temp[i,5:9],NA)
        if(!is.na(N.temp[i,11])){N.data.fit[indx,9]=N.temp[i,11]}
        indx=indx+1
    }
    if(N.temp[i,4]==0){
        N.data.fit[indx,] = c(N.temp[i,1:2],1,N.temp[i,5:9],NA)
        if(!is.na(N.temp[i,11])){N.data.fit[indx,9]=N.temp[i,11]}
        indx=indx+1
        N.data.fit[indx,] = c(N.temp[i,2:3],0,N.temp[i,5:9],NA)
        if(!is.na(N.temp[i,11])){N.data.fit[indx,9]=N.temp[i,11]}
        indx=indx+1
    }
}
N.data.fit
tail(N.data.fit,16)
#indexing records
n.fit=dim(N.data.fit)[1]

#filling in known CWD + status in column 9 of N.data.fit, which is equivalent to known z's in the analysis, remaining z's must be predicted
# for(i in 1:n.fit){
#     if(is.na(N.data.fit[i,9]) & N.data.fit[i,4]==1)N.data.fit[i,4]=1
# }

#for fast morts, remove the "living"contribution to survival lines
mort.check=which(N.data.fit[,3]==0)-1
rm.indx=c()
for(i in mort.check){
    if(N.data.fit[i,1]>=N.data.fit[i,2])rm.indx=c(rm.indx,i)
}

N.data.fit[sort(c(rm.indx,rm.indx+1)),]
rm.indx
N.data.fit=N.data.fit[-rm.indx,]

# indexing records
n.fit=dim(N.data.fit)[1]
n.fit

#replace those bow harvested labeled with gun harvest, as 0...

N.data.fit[which(N.data.fit[,8]==1),6]=0

#set initial values for z's
z = N.data.fit[,9]
z.indx = which(is.na(z))
n.z = length(z.indx)
# z.init=N.data.fit[,9]
# z.init[is.na(z.init)]=0
sum(N.data.fit[,4])
sum(N.data.fit[,9],na.rm=TRUE)

tail(N.data.fit,40)



######################################################################################################################3
###
### Adding covariates i.e. CWD status at capture and hunter harvest
###
######################################################################################################################3

###
### define Jags model 
###

sink("model.fawn.R")
cat("
    model{
    
    # Priors
    beta0 ~ dnorm(0, 1.0E-6) # intercept log hazard rate
    beta3 ~ dnorm(0, 1.0E-6) # gun harvest effect log hazard rate


    # Likelihood for the total hazard
    for (j in 1:records) {
        for (k in left[j]:(right[j]-1)) {
          UCH[j,k] <- exp(beta0 +beta3*x5[k] + timeuse[k]) #+beta2*x4[k]
        }

        # total prob of surviving
        SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
        censor[j] ~ dbern(SLR[j])      
    }

    #Auto-regressive random effect for time
    time[1] ~ dnorm(0,.001)
    for(t in 2:nT){
    time[t] ~ dnorm(time[t-1],tau)
    }
    #sum to zero constraint
    for(t in 1:nT){
      timeuse[t]=time[t] - mean(time)
    }
        
    tau <- sd^-2
    sd ~ dunif(0,10)

    #Derived parameters
    for (t in 1:nT){
        llambda.out[t]<-beta0 + beta3*x5[t] + timeuse[t] #CWD- gun harvest only
        UCH0[t]<-exp(llambda.out[t])
        CH0[t]<-sum(UCH0[1:t])
        S0[t]<-exp(-CH0[t])
    }
}
    ",fill = TRUE)
sink()

#specify initial values
inits1 <- list(beta0=-5,beta3=-1, 
               .RNG.name = "base::Wichmann-Hill", .RNG.seed=1)
inits2 <- list(beta0=-7,beta3=-1, 
               .RNG.name = "base::Marsaglia-Multicarry", .RNG.seed=2)
inits3 <- list(beta0=-8,beta3=-3, 
               .RNG.name = "base::Super-Duper", .RNG.seed=3)
inits <- list(inits1, inits2, inits3)


#identify params to monitor
parameters<-c("beta0","beta3","S0","sd","CH0") #"beta2","beta3",

###
### n.fit = # records, dim(N.data.fit)[1]
### x1 = cwd affect <- beta1 <- N.temp[,5].... N.data.fit[,4],"ps
### x4 = for derived parameters, bow season survival
### x5 = for derived parameters, gun season survival
### z = N.data.fit[,9], KNOWN cwd status from post mortem necropsy test
### z.indx = z's that need to be drawn
### n.z = # of z's that need to be drawn
###

# Bundle data
jags.data <- list(records = n.fit,
                  left = N.data.fit[,1],
                  right = N.data.fit[,2],
                  censor = N.data.fit[,3],
                  nT = study.end,
                  x5 = gun.harvest.haz
                  )

# x4 = bow.harvest.haz,

# MCMC settings
nt = 10
ni = 200000
nb = 200000
nc = 3

# call parallel version of jags using dclone
# cl <- makeCluster(3,type="SOCK")
# clusterEvalQ(cl,library(dclone))
# clusterEvalQ(cl,setwd(getwd()))
# out=jags.parfit(cl, data=jags.data, params=parameters, model="model.fawn.R",inits=inits,n.chains=nc,n.thin=10,n.iter=ni,n.burnin=nb)
# inits <- list(inits1, inits2, inits3)
# clusterExport(cl, c("jags.data", "params", "inits"))
# 
# jagsparallel <- function(i,...){
#     jags.fit(data=jags.data,params=params,model="model.fawn.R",inits=inits[[i]],
#              n.chains=1,updated.model=FALSE,...)
# }
# 
# out <- parLapply(cl,1:nc,jagsparallel,n.adapt=10000,n.update=nb,n.iter=ni,thin=nt)
# time <- Sys.time() - start
# stopCluster(cl)
# out <- as.mcmc.list(lapply(out,as.mcmc))
cl<-makeCluster(3,timeout=5184000)
clusterSetRNGStream(cl = cl, iseed = 121915)
clusterExport(cl, c("jags.data","parameters"))
for (i in seq_along(cl)){
    init <- inits[[i]]
    clusterExport(cl[i], "init")
}


out<-clusterEvalQ(cl, {
    library(rjags)
    load.module("dic")
    #load.module("glm")
    out1<-jags.model("model.fawn.R",jags.data,init,n.adapt=10000,n.chains=1)
    out2<-coda.samples(out1,parameters,n.iter=200000,thin=20)
    return(as.mcmc(out2))
})
class(out)
out <- as.mcmc(out)
traceplot(out[,"beta0"])
traceplot(out[,"beta3"])
traceplot(out[,"sd"])
gelman.diag(out)

#stopCluster(cl)

# If has not converged, continue updating:
# out<-clusterEvalQ(cl, {
#     out2<-coda.samples(out1,parameters,n.iter=200,thin=1)
#     return(as.mcmc(out2))
# })
# out.mcmc <- as.mcmc(out)
# traceplot(out.mcmc[,"sigma"])
# stopCluster(cl)

n.days=study.end
fit.sum=summary(out)[[1]]
fit.quant=summary(out)[[2]]
tail(fit.sum)
tail(fit.quant)

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
CH0.mn=fit.sum[1:n.days,1]
CH0.l=fit.quant[1:n.days,1]
CH0.u=fit.quant[1:n.days,5]

out.ch0= data.frame(Days,CH0.mn,CH0.l,CH0.u)

ggplot(data = out.ch0, aes(x = Days,y=CH0.mn))+
    geom_ribbon(aes(ymin=CH0.l,ymax=CH0.u),alpha=.1,show.legend=NA,linetype=0)+
    geom_line()+theme_bw()
     

fit.sum[n.days:(n.days+10),]
tail(fit.sum)     

Survival=fit.sum[(n.days+1):(2*n.days),1]
Lower=fit.quant[(n.days+1):(2*n.days),1]
Upper=fit.quant[(n.days+1):(2*n.days),5]

out.df=data.frame(cbind(Days,Survival,Lower,Upper))


ggplot(data = out.df, aes(x = Days,y=Survival))+
    geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=.1,show.legend=NA,linetype=0)+
    geom_line()+theme_bw()

###
### Hazard analysis coefficients
### 

output=cbind(fit.sum[(2*n.weeks+1):(2*n.weeks+2),1:2],fit.quant[(2*n.weeks+1):(2*n.weeks+2),c(1,5)])
row.names(output)=c("Intercept","Gun Harvest")
output=round(output,2)
output

write.csv(output,file="Hazard_table_fawn.csv",row.names = TRUE)

print.xtable(xtable(output,caption = 'Coefficients of the hazard rate analysis where the baseline hazard shows the daily hazard intercept, and additive term for Gun Harvests. The first column denotes the mean of the posterior distribution with the standard deviation (SD), and equal-tailed Bayesian credible intervals in the right two columns.',
             align = c('c','c','c','c','c'),label="tab:hazard_fawn",digits=c(0,2,2,2,2)),
             sanitize.text.function = function(x) {x},
             hline.after=c(0,2),
             include.rownames=TRUE,
             rotate.colnames = FALSE,
             caption.placement = "top")

######################################################################################################################3
###
### Summary Output
###
######################################################################################################################3


###
### N.temp[,1] = e
### N.temp[,2] = r 
### N.temp[,3] = s
### N.temp[,4] = censor
### N.temp[,5] = CWD status
### N.temp[,6] = body condition
### N.temp[,7] = gun harvest
### N.temp[,8] = sex
### N.temp[,9] = bow harvest
###

###################################################################################################
###
### Harvest summary plots
###
###################################################################################################

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
gun1.end=gun.end.week
gun2.end=gun.holiday.end.week

fit.sum[c(study.end+bow.end),1]

post.bow.harvest=cbind(fit.sum[study.end+bow.end,1],fit.sum[study.end+bow.end,2],fit.quant[study.end+bow.end,1],fit.quant[study.end+bow.end,5])
post.bow.harvest=round(post.bow.harvest,3)
post.bow.harvest
post.bow.harvest=data.frame(post.bow.harvest)
names(post.bow.harvest)=c("Mean","SD","0.025","0.975")


post.gun1.harvest=cbind(fit.sum[study.end+gun1.end,1],fit.sum[study.end+gun1.end,2],fit.quant[n.weeks+gun1.end,1],fit.quant[study.end+gun1.end,5])
post.gun1.harvest=round(post.gun1.harvest,3)
post.gun1.harvest
post.gun1.harvest=data.frame(post.gun1.harvest)
names(post.gun1.harvest)=c("Mean","SD","0.025","0.975")


post.gun2.harvest=cbind(fit.sum[study.end+gun2.end,1],fit.sum[study.end+gun2.end,2],fit.quant[study.end+gun2.end,1],fit.quant[study.end+gun2.end,5])
post.gun2.harvest=round(post.gun2.harvest,3)
post.gun2.harvest
post.gun2.harvest=data.frame(post.gun2.harvest)
names(post.gun2.harvest)=c("Mean","SD","0.025","0.975")



Period=c("Pre-Gun Harvest","","Gun","","Holiday Gun","")

survival_all_sum=data.frame(cbind(Period,rbind(post.bow.harvest,post.gun1.harvest,post.gun2.harvest)))
survival_all_sum
names(survival_all_sum)=c("Period","Mean","SD","0.025","0.975")

write.csv(survival_all_sum,"survival_fawn.csv",row.names=FALSE)

print.xtable(xtable(survival_all_sum,caption = 'The mean of the posterior distribution (black line) of the daily cumulative probability of survival of fawns after birth during 2017 of white-tailed deer in southwest Wisconsin, along with 95\\% equal-tailed Bayesian credible intervals (gray shaded region).', 
                    align = c('c','c','c','c','c','c'),label="tab:survival_fawn",digits=c(0,0,2,2,2,2)),
             sanitize.text.function = function(x) {x},
             hline.after=c(-1,-1,0,6,6),
             add.to.row = list(pos=list(c(2,4,6)), command="\\hdashline \n"),
             include.rownames=FALSE,
             rotate.colnames = FALSE,
             caption.placement = "top")  


yday('2017-12-31')

#day cut-off for harvest 
non.harvest.survival.end=yday('2017-09-16')-start.day
bow.start.day=yday('2017-09-16')-start.day+1
bow.end.day = yday('2018-01-07')+365-start.day+1
gun.start.day=yday('2017-11-18')-start.day+1
gun.end.day = yday('2017-11-26')-start.day+1
gun.holiday.start.day = yday('2017-12-24')-start.day+1
gun.holiday.end.day = yday('2018-01-01')+365-start.day+1
study.end=yday('2018-01-07')+365-start.day+1


pdf("Survival_fawn_v2.pdf", width=8, height =6)
ggplot(data =out.df,aes(x = Days,y=Survival))+
    annotate('rect',xmax=gun.holiday.end.day, xmin=gun.holiday.start.day, ymax=Inf, ymin=-Inf, linetype=2,alpha=.5,fill="grey90")+
    annotate('rect',xmax=gun.end.day, xmin=gun.start.day, ymax=Inf, ymin=-Inf, linetype=2,alpha=.5,fill="grey90")+
    geom_line()+theme_bw()+
    geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=.1,show.legend=NA,linetype=0)+
    ggtitle("Fawn Survival")+xlab("Time (Days)")+ylab("Survival Probability")+
    geom_vline(aes(xintercept=gun.start.day),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.end.day),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.holiday.start.day),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.holiday.end.day),linetype=2,color="grey50")+
    scale_x_discrete(limit = c(1,51,101,151,gun.start.day,gun.holiday.start.day),labels = c("May 16","July 5","Aug 24","Oct 13","Nov 18","Dec 24","Jan 7"))
dev.off()

# as.Date(start.day-1+50,origin="2017-01-01")#July 5
# as.Date(start.day-1+100,origin="2017-01-01")#Aug 24
# as.Date(start.day-1+150,origin="2017-01-01")#Oct 13


#saveworking directory
save.image("fawn_survival_v2.Rdata")

