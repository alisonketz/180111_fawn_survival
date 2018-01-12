
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
    
