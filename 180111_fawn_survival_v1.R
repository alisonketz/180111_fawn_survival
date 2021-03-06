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

load('fawn_survival_v7.Rdata')

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


#start of the study
start.week=week(ymd("2017-05-16"))

 
###
### removing Adults from the mortality data frame
###

ad.mort.rm=c()
for(l in 1:n.fawncap){
    ad.mort.rm=c(ad.mort.rm,which(d.mort$lowtag==d.cap$lowtag[l]))
}
d.fawnmort=d.mort[-ad.mort.rm,]

n.fawnmort=dim(d.fawnmort)[1]

###
### initial format of data - rows = individuals, cols = (e_i,r_i,s_i,censor)
###

#making dates discrete for fitting model

#e
d.cap$dateDis = (week(d.cap$date) - start.week + 1)
d.fawncap$dateDis = (week(d.fawncap$date) - start.week + 1)

#r
d.fawnmort$mORTAlertdateDis=(week(d.fawnmort$mORTAlertdate) - start.week + 1)
d.fawnmort$estMortDateDis = (week(d.fawnmort$estMortDate) - start.week + 1)

#s
d.fawnmort$collarFoundDis=(week(d.fawnmort$collarFound) - start.week + 1)



#week cut-off for harvest 
non.harvest.survival.end=week(ymd('2017-09-16'))-start.week
bow.start.week=week(ymd('2017-09-16'))-start.week+1
bow.end.week = week(ymd('2018-01-07'))+52-start.week+1
gun.start.week=week(ymd('2017-11-18'))-start.week+1
gun.end.week = week(ymd('2017-11-26'))-start.week+1
gun.holiday.start.week = week(ymd('2017-12-24'))-start.week+1
gun.holiday.end.week = week(ymd('2018-01-01'))+52-start.week+1
study.end=week(ymd('2018-01-07'))+52-start.week+1

d.fawnmort$bow=0
d.fawnmort$bow[which(d.fawnmort$weapon!="Rifle" & d.fawnmort$cause1 == "Hunter harvest")[-1]]=1
         
# likely gun harvest
which(d.fawnmort$weapon!="Rifle" & d.fawnmort$cause1 == "Hunter harvest")

#setting harvest.season for derived parameters
n.weeks=study.end

# non.harvest.survival.end
# bow.start.week
# bow.end.week
# gun.start.week
# gun.end.week
# gun.holiday.start.week
# gun.holiday.end.week


bow.harvest.hazard=bow.harvest.haz = c(rep(0,non.harvest.survival.end),rep(1,bow.end.week-bow.start.week+1))
bow.harvest.hazard
length(bow.harvest.haz)

gun.harvest.haz=c(rep(0,gun.start.week-1),rep(1,gun.end.week-gun.start.week+1),rep(0,gun.holiday.start.week-gun.end.week),rep(1,gun.holiday.end.week-gun.holiday.start.week))
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
    N.temp[i,1] = d.fawncap$dateDis[i]#e
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

#censor morts caused by capture
mort.killcap=which(d.fawnmort$cause1=="Capture Related")
mort.killcap
d.fawnmort[mort.killcap,]
cens.id=d.fawnmort$lowtag[mort.killcap]

for(i in 1:n.temp){
    for(j in cens.id){
        if(N.temp[i,10]==j){
            N.temp[i,2:4]=c(study.end,study.end,1)
        }
    }
}

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




#set harvest coefficient

# which(N.temp[,3]>non.harvest.survival.end)

# N.temp[which(N.temp[,3]>non.harvest.survival.end),7]=1
# which(N.temp[,3]>non.harvest.survival.end & N.temp[,4]==0)
# sum(N.temp[,7])#number of gun harvested
# N.temp[,7]
# 
# head(N.temp)
# N.temp[6,1:3]=c(1,2,2)

#making last known alive be week before mort
for(i in 1:n.temp){
    if(N.temp[i,4]==0 & N.temp[i,2]==N.temp[i,3]){
        N.temp[i,2]=N.temp[i,2]-1
    }
}

N.temp[which(N.temp[,2]<N.temp[,1])+1,1]=N.temp[which(N.temp[,2]<N.temp[,1])+1,2]




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

N.data.fit = matrix(NA,nr=n.temp+n.fawnmort-16,ncol=9)

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
    beta2 ~ dnorm(0, 1.0E-6) # bow harvest effect log hazard rate
    beta3 ~ dnorm(0, 1.0E-6) # gun harvest effect log hazard rate


    # Likelihood for the total hazard
    for (j in 1:records) {
        for (k in left[j]:(right[j]-1)) {
            UCH[j,k] <- exp(beta0 +beta2*x4[k]+beta3*x5[k])
        }
        
        # total prob of survivingy
        SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
        censor[j] ~ dbern(SLR[j])      
    }

    #Derived parameters
    for (t in 1:n.weeks){
        llambda.out[t]<-beta0 + beta2*x4[t] +beta3*x5[t]#CWD-
        UCH0[t]<-exp(llambda.out[t])
        CH0[t]<-sum(UCH0[1:t])
        S0[t]<-exp(-CH0[t])
    }

}
    ",fill = TRUE)
sink()

#specify initial values
inits<-list(list("beta0"=1,"beta2"=1,"beta3"=.1),
            list("beta0"=-.5,"beta2"=2,"beta3"=-.3),
            list("beta0"=1,"beta2"=.5,"beta3"=1.6)
            )

#identify params to monitor
parameters<-c("beta0","beta2","beta3","S0") 

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
                  n.weeks = study.end,
                  x4 = bow.harvest.haz,
                  x5 = gun.harvest.haz
                  )


# MCMC settings
nt = 1
ni = 10000
nb = 10000
nc = 3

# call parallel version of jags using dclone
cl=makeCluster(3)
out.cwd=jags.parfit(cl, data=jags.data, params=parameters, model="model.fawn.R",inits=inits,n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
stopCluster(cl)

# gelman.diag(out.cwd)
n.weeks=max(N.data.fit,na.rm=TRUE)

fit.sum=summary(out.cwd)[[1]]
fit.quant=summary(out.cwd)[[2]]
fit.sum
fit.quant


n.sum=dim(fit.sum)[1]
round(fit.sum[(n.sum-188):n.sum,1])
which(N.data.fit[,4]==1)
sum(N.data.fit[,4])
z.out=as.numeric(fit.sum[(n.sum-188):n.sum,1]>.2)
sum(z.out)

which(z.out==1)


# The Colorblind palette with grey:
cbPalette <- c( "#0072B2", "#D55E00", "#CC79A7","#999999", "#E69F00", "#56B4E9")

###
### plotting single survival curve
###

Weeks=rep(1:n.weeks)
Survival=fit.sum[1:(n.weeks),1]
Lower=fit.quant[1:(n.weeks),1]
Upper=fit.quant[1:(n.weeks),5]

out=data.frame(cbind(Weeks,Survival,Lower,Upper))



Weeks=rep(1:n.weeks,2)
CWD.status=(c(rep(0,n.weeks),rep(1,n.weeks)))

Survival=fit.sum[1:(2*n.weeks),1]
Lower=fit.quant[1:(2*n.weeks),1]
Upper=fit.quant[1:(2*n.weeks),5]

out=data.frame(cbind(Weeks,CWD.status,Survival,Lower,Upper))
out$CWD.status=as.factor(CWD.status)


pdf("Fawn_survival_v1.pdf",width=8,height=6)
ggplot(data = out, aes(x = Weeks,y=Survival))+
    annotate('rect',xmax=gun.holiday.end.week, xmin=gun.holiday.start.week, ymax=Inf, ymin=-Inf, linetype=2,alpha=.5,fill="grey90")+
    annotate('rect',xmax=gun.end.week, xmin=gun.start.week-1, ymax=Inf, ymin=-Inf, linetype=2,alpha=.5,fill="grey90")+
    geom_line()+theme_bw()+
    geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=.1,show.legend=NA,linetype=0)+
    ggtitle("Fawn Survival")+xlab("Time")+ylab("Survival Probability")+
    geom_vline(aes(xintercept=bow.start.week-1),linetype=2,color="grey70")+
    geom_vline(aes(xintercept=gun.start.week-1),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.end.week),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.holiday.start.week),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.holiday.end.week),linetype=2,color="grey50")+
    geom_text(x=36.85,y=1,label="Bow",color="grey70")+
    geom_text(x=45.5,y=1,label="Gun",color="grey50")+
    scale_x_discrete(limit = c(1,9,17,26,34),labels = c("May","July","Sep 16","Nov 18","Jan 7"))
dev.off()


bow.end.week
###
### Hazard analysis coefficients
### 

output=cbind(fit.sum[(2*n.weeks+1):(2*n.weeks+4),1:2],fit.quant[(2*n.weeks+1):(2*n.weeks+4),c(1,5)])
row.names(output)=c("CWD(-), No harvest","CWD(+), No harvest","Bow Harvest","Gun Harvest")
output=round(output,2)
output
write.csv(output,file="Hazard_table_v7.csv",row.names = TRUE)

print.xtable(xtable(output,caption = 'Coefficients of the hazard rate analysis where the baseline hazard shows the hazard for chronic wasting disease positive adult individuals not exposed to harvest, with additive hazard effects in the following three rows. The first column denotes the mean of the posterior distribution with the standard deviation (SD), and equal-tailed Bayesian credible intervals in the right two columns.', 
                    align = c('c','c','c','c','c'),label="tab:hazard_v7",digits=c(0,0,3,3,3)),
             sanitize.text.function = function(x) {x},
             hline.after=c(0,4),
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



###
### Frequency Tables, EDA
###


eda.out=data.frame(N.temp)
names(eda.out)=c("Enter","AliveWeek","Exit","Censor","CWDstatus","BodyCondition","Gun","Sex","Bow")

pos.eda.out=eda.out[eda.out$CWDstatus==1,]
neg.eda.out=eda.out[eda.out$CWDstatus==0,]

apply(neg.eda.out,2,sum)[c(4,7,8,9)]
apply(pos.eda.out,2,sum)[c(4,7,8,9)]


freq.table=data.frame(matrix(NA,nr=2,ncol=5))
freq.table[,1]=c("Negative","Positive")
names(freq.table)=c("CWD.status","Alive","NumberFemale","GunHarvest","BowHarvest")
freq.table[1,2:5]=apply(neg.eda.out,2,sum)[c(4,8,7,9)]
freq.table[2,2:5]=apply(pos.eda.out,2,sum)[c(4,8,7,9)]



myfreq=data.frame(t(as.matrix(freq.table[,2:5])))
names(myfreq)=c("Negative","Positive")
pdf("CWD_freq_table.pdf")
barplot(as.matrix(myfreq), main="CWD Status Frequency Table", ylab="Total", beside=TRUE, col=cbPalette[3:6])
legend(6, 80, c("Number Collars Alive","Number Collared Female","Gun Harvested","Bow Harvested"), fill=cbPalette[3:6])
dev.off()

###
### No "Alive" Category
###

freq.table=data.frame(matrix(NA,nr=2,ncol=4))
freq.table[,1]=c("Negative","Positive")
names(freq.table)=c("CWD.status","NumberFemale","GunHarvest","BowHarvest")
freq.table[1,2:4]=apply(neg.eda.out,2,sum)[c(8,7,9)]
freq.table[2,2:4]=apply(pos.eda.out,2,sum)[c(8,7,9)]

myfreq=data.frame(t(as.matrix(freq.table[,2:4])))
names(myfreq)=c("Negative","Positive")
pdf("CWD_freq_table_v2.pdf")
barplot(as.matrix(myfreq), main="CWD Status Frequency Table", ylab="Total", beside=TRUE, col=cbPalette[4:6])
legend(5, 50, c("Number Collared Female","Gun Harvested","Bow Harvested"), fill=cbPalette[4:6])
dev.off()

###
### Males vs females, breakdown
###

male.eda.out=eda.out[eda.out$Sex==1,]
female.eda.out=eda.out[eda.out$Sex==0,]


freq.table=data.frame(matrix(NA,nr=2,ncol=5))
freq.table[,1]=c("Male","Female")
names(freq.table)=c("Sex","Alive","CWD.status","GunHarvest","BowHarvest")
freq.table[1,2:5]=apply(male.eda.out,2,sum)[c(4,5,7,9)]
freq.table[2,2:5]=apply(female.eda.out,2,sum)[c(4,5,7,9)]


myfreq=data.frame(t(as.matrix(freq.table[,2:5])))
myfreq=myfreq[,2:1]
names(myfreq)=c("Females","Males")
myfreq
pdf("Sex_freq_table.pdf")
barplot(as.matrix(myfreq), main="Sex Frequency Table", ylab="Total", beside=TRUE, col=cbPalette[3:6])
legend(7, 62, c("Alive","CWD(+) at Capture", "Gun Harvested","Bow Harvested"), fill=cbPalette[3:6])
dev.off()











###################################################################################################
###
### Harvest summary plots
###
###################################################################################################


N.harvest=N.temp[which(N.temp[,7]==1),]
N.harvest[which( N.harvest[,9]==1),7]=0

total.harvest=dim(N.harvest)[1]

male.harvest=N.harvest[N.harvest[,8]==1,]
female.harvest=N.harvest[N.harvest[,8]==0,]

dim(N.harvest)[1]-sum(N.harvest[,8])

male.harvest[,5]

harvest.tab=data.frame(matrix(NA,nr=2,ncol=5))

harvest.tab[,1]=c("Male","Female")
names(harvest.tab)=c("Sex","TotalHarvest","CWD.status","GunHarvest","BowHarvest")
harvest.tab[1,2:5]=apply(male.harvest,2,sum)[c(8,5,7,9)]
harvest.tab[2,2:5]=apply(female.harvest,2,sum)[c(8,5,7,9)]
harvest.tab[2,2]=dim(female.harvest)[1]

harvest.tab


myfreq=data.frame(t(as.matrix(harvest.tab[,2:5])))
myfreq=myfreq[,2:1]
names(myfreq)=c("Females","Males")
myfreq
pdf("Harvest_freq_table.pdf")
barplot(as.matrix(myfreq), main="Harvested Collars Frequency Table", ylab="Total", beside=TRUE, col=cbPalette[3:6])
legend(1.5, 10, c("Total Harvest","CWD(+) at Capture", "Gun Harvested","Bow Harvested"), fill=cbPalette[3:6])
dev.off()





###################################################################################################
###
### the "annual" survival estimates (i.e., the mean/median and 95%CI for surviving to the end of the analysis period) broken down by harvest types
###
###################################################################################################

###
### Table with survival before/after bow and gun hunts
###
pre.end=non.harvest.survival.end
bow.end=gun.start.week-1
gun1.end=gun.end.week
gun2.end=gun.holiday.end.week

fit.sum[c(pre.end,n.weeks+pre.end),1]

pre.harvest=cbind(fit.sum[c(pre.end,n.weeks+pre.end),1],fit.sum[c(pre.end,n.weeks+pre.end),2],fit.quant[c(pre.end,n.weeks+pre.end),1],fit.quant[c(pre.end,n.weeks+pre.end),5])
pre.harvest=round(pre.harvest,3)
pre.harvest
pre.harvest=data.frame(pre.harvest)
pre.harvest=cbind(rep(NA,2),pre.harvest)
pre.harvest[,1]=c("Negative","Positive")
names(pre.harvest)=c("CWD Status","Mean","SD","0.025","0.975")



post.bow.harvest=cbind(fit.sum[c(bow.end,n.weeks+bow.end),1],fit.sum[c(bow.end,n.weeks+bow.end),2],fit.quant[c(bow.end,n.weeks+bow.end),1],fit.quant[c(bow.end,n.weeks+bow.end),5])
post.bow.harvest=round(post.bow.harvest,3)
post.bow.harvest
post.bow.harvest=data.frame(post.bow.harvest)
post.bow.harvest=cbind(rep(NA,2),post.bow.harvest)
post.bow.harvest[,1]=c("Negative","Positive")
names(post.bow.harvest)=c("CWD Status","Mean","SD","0.025","0.975")


post.gun1.harvest=cbind(fit.sum[c(gun1.end,n.weeks+gun1.end),1],fit.sum[c(gun1.end,n.weeks+gun1.end),2],fit.quant[c(gun1.end,n.weeks+gun1.end),1],fit.quant[c(gun1.end,n.weeks+gun1.end),5])
post.gun1.harvest=round(post.gun1.harvest,3)
post.gun1.harvest
post.gun1.harvest=data.frame(post.gun1.harvest)
post.gun1.harvest=cbind(rep(NA,2),post.gun1.harvest)
post.gun1.harvest[,1]=c("Negative","Positive")
names(post.gun1.harvest)=c("CWD Status","Mean","SD","0.025","0.975")



# 
# post.bow2.harvest=cbind(fit.sum[c(bow2.end,n.weeks+bow2.end),1],fit.sum[c(bow2.end,n.weeks+bow2.end),2],fit.quant[c(bow2.end,n.weeks+bow2.end),1],fit.quant[c(bow2.end,n.weeks+bow2.end),5])
# post.bow2.harvest=round(post.bow2.harvest,3)
# post.bow2.harvest
# post.bow2.harvest=data.frame(post.bow2.harvest)
# post.bow2.harvest=cbind(rep(NA,2),post.bow2.harvest)
# post.bow2.harvest[,1]=c("Negative","Positive")
# names(post.bow2.harvest)=c("CWD Status","Mean","SD","0.025","0.975")


post.gun2.harvest=cbind(fit.sum[c(gun2.end,n.weeks+gun2.end),1],fit.sum[c(gun2.end,n.weeks+gun2.end),2],fit.quant[c(gun2.end,n.weeks+gun2.end),1],fit.quant[c(gun2.end,n.weeks+gun2.end),5])
post.gun2.harvest=round(post.gun2.harvest,3)
post.gun2.harvest
post.gun2.harvest=data.frame(post.gun2.harvest)
post.gun2.harvest=cbind(rep(NA,2),post.gun2.harvest)
post.gun2.harvest[,1]=c("Negative","Positive")
names(post.gun2.harvest)=c("CWD Status","Mean","SD","0.025","0.975")



Period=c("Pre-Harvest","","Bow/Archery","","Gun","","Holiday Gun","")

survival_all_sum=data.frame(cbind(Period,rbind(pre.harvest,post.bow.harvest,post.gun1.harvest,post.gun2.harvest)))
survival_all_sum
names(survival_all_sum)=c("Period","CWD Status","Mean","SD","0.025","0.975")

write.csv(survival_all_sum,"survival_all_sum_v7.csv",row.names=FALSE)

print.xtable(xtable(survival_all_sum,caption = 'survival probability nlahshfls', 
                    align = c('c','c','c','c','c','c','c'),label="tab:survival_v7",digits=c(0,0,2,2,2,2,2)),
             sanitize.text.function = function(x) {x},
             hline.after=c(-1,-1,0,8,8),
             add.to.row = list(pos=list(c(2,4,6)), command="\\hdashline \n"),
             include.rownames=FALSE,
             rotate.colnames = FALSE,
             caption.placement = "top")  

###
### Prevalence posterior plot
###

out.psi.indx=dim(fit.sum)[1]
fit.df=data.frame(rbind(out.cwd[[1]],out.cwd[[2]],out.cwd[[3]])[,out.psi.indx])

plot(density(fit.df$Psi))
names(fit.df)="Psi"

cols=c("grey50","darkred")

pdf("Prevalence.pdf")
ggplot(data=fit.df,aes(x=Psi))+geom_density(size=1)+
    geom_vline(aes(xintercept=fit.sum[out.psi.indx,1],color=cols[2]),linetype=2)+
    geom_vline(aes(xintercept=fit.quant[out.psi.indx,1],color=cols[1]),linetype=2)+
    geom_vline(aes(xintercept=fit.quant[out.psi.indx,5],color=cols[1]),linetype=2)+
    scale_color_manual(values=cols,name="Posterior Statistics",labels=c("Mean = .15","Credible Interval = (.10,.21)"))+
    ggtitle("Prevalence")+xlab(expression(psi))+ylab("Density")+
    theme_bw()
dev.off()

prevalence.out=c(fit.sum[out.psi.indx,1],fit.quant[out.psi.indx,1],fit.quant[out.psi.indx,5])





#saveworking directory
save.image("survival_v7.Rdata")

