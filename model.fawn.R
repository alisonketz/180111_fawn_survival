
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
    
