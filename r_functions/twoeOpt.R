twoeOpt<-function(model, data){
  cost<-function(x){
      yhat<-data.frame(Pred=numeric(), Pcorr2=numeric())
      for(i in unique(data$C_AP)){
          out<-as.data.frame(ode(y=c(P=0, S=i-mean(data[(data$C_AP==i), "Pblanks"], na.rm=T)), 
                                 parms = c(Vmax1=x[1], 
                                           Vmax2=x[2],
                                           Km=x[3], Kic=x[4]), 
                                 model, times=unique(data[data$C_AP==i, "Time"]))) 
          out<-out[, c("time", "P")]
          colnames(out)<-c("Time", "Pred")
          outm<-merge(out, data[(data$C_AP==i), c("Time", "Pcorr2")], by = c("Time"))[,-1]
          yhat<-rbind(yhat, outm)
        
        
      }
      
      RMSE<-with(yhat, sum((((Pred-Pcorr2))^2), na.rm = T))
      return(RMSE)
    }
    
    par_mcmc<-modMCMC(f=cost, p=c(1e-4, 1e-2, 50, 30), 
                      lower=c(1e-6, 1e-6, 1e-6, 1e-6),
                      upper=c(1e3, 1e3, 1e3, 1e3), niter=10000)
    #lower and upper limits for parameters are extracted
    pl<-summary(par_mcmc)["min",]
    pu<-summary(par_mcmc)["max",]
    
    #these limits are used to find global optimum by DEoptim
    opt_par<-DEoptim(fn=cost, lower=pl, upper=pu, 
                     control = c(itermax = 10000, steptol = 50, reltol = 1e-8, 
                                 trace=FALSE, strategy=3, NP=250))
    rsq<-function(x){
      yhat<-data.frame(Time = numeric(), Pred=numeric(), Pcorr2=numeric(), C_AP=numeric())
      for(i in unique(data$C_AP)){
        out<-as.data.frame(ode(y=c(P=0, S=i-mean(data[(data$C_AP==i), "Pblanks"], na.rm=T)), 
                               parms = c(Vmax1=x[1], 
                                         Vmax2=x[2],
                                         Km=x[3], Kic=x[4]), 
                               model, times=unique(data[data$C_AP==i, "Time"])))
          out<-out[, c("time", "P")]
          colnames(out)<-c("Time", "Pred")
          outm<-merge(out, data[(data$C_AP==i), c("Time", "Pcorr2")], by = c("Time"))
          outm$C_AP<-rep(mean(as.numeric(data[(data$C_AP==i), "C_AP"])), times = nrow(outm))
          yhat<-rbind(yhat, outm)
        
        
      }
      
      SSres=with(yhat, sum(((Pcorr2-Pred)^2), na.rm = T))
      SStot=with(yhat, sum(((Pcorr2-mean(Pcorr2, na.rm = T))^2), na.rm = T))
      ll=with(yhat, -sum(((Pcorr2-Pred)^2), na.rm = T)/2/(sd(Pcorr2, na.rm = T)^2))
      R2<-1-SSres/SStot
      N<-length(x)
      AIC<-2*N-2*ll
      Gfit<-c(R2=R2, N=N, AIC=AIC, ll=ll, SSres=SSres, SStot=SStot)
      rsq_out<-list(Yhat=yhat, Gfit=Gfit)
      return(rsq_out)
    }
    
    out_all<-list(P = opt_par$optim$bestmem,
                  Goodness = rsq(as.numeric(opt_par$optim$bestmem)),
                  MCMC = par_mcmc)
    return(out_all)
    
}