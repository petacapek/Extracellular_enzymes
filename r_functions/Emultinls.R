Emultinls<-function(data){
  #1. Define functions to evaluate
  ###Without inhibition
  wi<-"Time ~ -1/Vmax*(Km*log((C_AP-Pblanks-Pcorr2)/(C_AP-Pblanks))-Pcorr2)"
  wic<-function(x){
    return(sum((d$Time - -1/x[1]*(x[2]*log((d$C_AP-d$Pblanks-d$Pcorr2)/(d$C_AP-d$Pblanks))-d$Pcorr2))^2))
  }
  ###Substrate inhibition
  si<-"Time ~ (Km*log((C_AP-Pblanks)/(C_AP-Pblanks-Pcorr2))+Pcorr2+((C_AP-Pblanks)^2-(C_AP-Pblanks-Pcorr2)^2)/2/Ks)/Vmax"
  sic<-function(x){
    return(sum((d$Time - (x[2]*log((d$C_AP-d$Pblanks)/(d$C_AP-d$Pblanks-d$Pcorr2))+d$Pcorr2+((d$C_AP-d$Pblanks)^2-(d$C_AP-d$Pblanks-d$Pcorr2)^2)/2/x[3])/x[1])^2))
  }
  ###Non-competitive inhibition
  nci<-"Time ~ -1/Vmax*(Km*((C_AP-Pblanks)/Ki+1)*log((C_AP-Pblanks-Pcorr2)/(C_AP-Pblanks))+
                      (1+Km/Ki+(C_AP-Pblanks)/Ki)*-Pcorr2-((C_AP-Pblanks-Pcorr2)^2-(C_AP-Pblanks)^2)/2/Ki)"
  ncic<-function(x){
    return(sum((d$Time - -1/x[1]*(x[2]*((d$C_AP-d$Pblanks)/x[3]+1)*log((d$C_AP-d$Pblanks-d$Pcorr2)/(d$C_AP-d$Pblanks))+
                                    (1+x[2]/x[3]+(d$C_AP-d$Pblanks)/x[3])*-d$Pcorr2-((d$C_AP-d$Pblanks-d$Pcorr2)^2-(d$C_AP-d$Pblanks)^2)/2/x[3]))^2))
  }
  ###Competitive inhibition
  ##b) competitive inhibition
  ci<-"Time ~ -1/Vmax*(Km*((C_AP-Pblanks)/Ki+1)*log((C_AP-Pblanks-Pcorr2)/(C_AP-Pblanks))+
                      (1+Km/Ki)*-Pcorr2)"
  cic<-function(x){
    return(sum((d$Time- -1/x[1]*(x[2]*((d$C_AP-d$Pblanks)/x[3]+1)*log((d$C_AP-d$Pblanks-d$Pcorr2)/(d$C_AP-d$Pblanks))+
                                   (1+x[2]/x[3])*-d$Pcorr2))^2))
  }
  ###Uncompetitive inhibition
  uci<-"Time ~ -1/Vmax*(Km*log((C_AP-Pblanks-Pcorr2)/(C_AP-Pblanks))+
                      (1+Km/Ki+(C_AP-Pblanks)/Ki)*-Pcorr2-((C_AP-Pcorr2-Pblanks)^2-(C_AP-Pblanks)^2)/2/Ki)"
  ucic<-function(x){
    return(sum((d$Time - -1/x[1]*(x[2]*log((d$C_AP-d$Pblanks-d$Pcorr2)/(d$C_AP-d$Pblanks))+
                                    (1+x[2]/x[3]+(d$C_AP-d$Pblanks)/x[3])*-d$Pcorr2-((d$C_AP-d$Pcorr2-d$Pblanks)^2-(d$C_AP-d$Pblanks)^2)/2/x[3]))^2))
  }
  #2. initial parameters
  swi<-data.frame(Vmax=c(1e-3, 1e2), Km=c(0.1, 500))
  ssi<-data.frame(Vmax=c(1e-3, 1e2), Km=c(0.1, 500), Ks=c(0.1, 500))
  snci<-data.frame(Vmax=c(1e-3, 1e2), Km=c(0.1, 500), Ki=c(1e-3, 500))
  swil<-c(1e-3, 0.1)
  swiu<-c(1e2, 500)
  ssil<-c(1e-3, 0.1, 1e-3)
  ssiu<-c(1e2, 500, 500)
  sncil<-c(1e-3, 0.1, 1e-3)
  snciu<-c(1e2, 500, 500)
  
  #3. For each slurry separatelly - run
  ##Initialize data frame where results will be stored
  sep_out<-data.frame(Model=character(), Dilution=numeric(), Vmax=numeric(), Km=numeric(),
                      Ks=numeric(), Ki=numeric(),
                      SSres=numeric(), SStot=numeric(), P=numeric(),
                      n=numeric())
  
  ##run
  for(i in unique(data$D)){
    ###Without inhibition
    tryCatch({
      wia<-nls2(wi, data=data[data$D==i, ], start=swi, control=list(maxiter=1e5), algorithm="brute-force")
      wiout<-nls2(wi, data=data[data$D==i, ], start=wia, control=list(maxiter=1e5, tol=1e-2, minFactor=1e-12))
      sep_out<-rbind(sep_out,
                     data.frame(Model=c("WI"),
                                Dilution=i, 
                                Vmax=coef(wiout)[1], 
                                Km=coef(wiout)[2],
                                Ks=NA, Ki=NA, 
                                SSres=sum(resid(wiout)^2), 
                                SStot=with(data[data$D==i, ], sum((Time-mean(Time))^2)), 
                                P=2, n=nrow(data[data$D==i, ])))
    }, error = function(e){print("NLS cannot converge - trying Differential Optimization algorithm")},
    finally = {
      d = data[data$D==i, ]
      wiout<-DEoptim(lower=swil, upper=swiu, fn=wic,
                      control = c(itermax = 10000, 
                                  trace=FALSE, strategy=3, NP=250))
      sep_out<-rbind(sep_out,
                     data.frame(Model=c("WI"),
                                Dilution=i, 
                                Vmax=as.numeric(wiout$optim$bestmem[1]), 
                                Km=as.numeric(wiout$optim$bestmem[2]),
                                Ks=NA, Ki=NA,
                                SSres=as.numeric(wiout$optim$bestval), 
                                SStot=with(data[data$D==i, ], sum((Time-mean(Time))^2)), 
                                P=2, n=nrow(data[data$D==i, ])))
    })
    
    ###Substrate inhibition
    tryCatch({
      sia<-nls2(si, data=data[data$D==i, ], start=ssi, control=list(maxiter=1e5), algorithm="brute-force")
      si_out<-nls2(si, data=data[data$D==i, ], start=sia, control=list(maxiter=1e5, tol=1e-2, minFactor=1e-12))
      sep_out<-rbind(sep_out,
                     data.frame(Model=c("SI"),
                                Dilution=i, 
                                Vmax=coef(si_out)[1], 
                                Km=coef(si_out)[2],
                                Ks=coef(si_out)[3], 
                                Ki=NA, 
                                SSres=sum((resid(si_out)^2)), 
                                SStot=with(data[data$D==i, ], sum((Time-mean(Time))^2)), 
                                P=3, n=nrow(data[data$D==i, ])))
    }, error = function(e){print("NLS cannot converge - trying Differential Optimization algorithm")},
    finally = {
      d = data[data$D==i, ]
      si_out<-DEoptim(lower=ssil, upper=ssiu, fn=sic,
                      control = c(itermax = 10000, 
                                  trace=FALSE, strategy=3, NP=250))
      sep_out<-rbind(sep_out,
                     data.frame(Model=c("SI"),
                                Dilution=i, 
                                Vmax=as.numeric(si_out$optim$bestmem[1]), 
                                Km=as.numeric(si_out$optim$bestmem[2]),
                                Ks=as.numeric(si_out$optim$bestmem[3]), 
                                Ki=NA, 
                                SSres=as.numeric(si_out$optim$bestval), 
                                SStot=with(data[data$D==i, ], sum((Time-mean(Time))^2)), 
                                P=3, n=nrow(data[data$D==i, ])))
    })
    
    ###Non-competitive inhibition
    tryCatch({
      ncia<-nls2(nci, data=data[data$D==i, ], start=snci, control=list(maxiter=1e5), algorithm="brute-force")
      nci_out<-nls2(nci, data=data[data$D==i, ], start=ncia, control=list(maxiter=1e5, tol=1e-2, minFactor=1e-12))
      sep_out<-rbind(sep_out,
                     data.frame(Model=c("NCI"),
                                Dilution=i, 
                                Vmax=coef(nci_out)[1], 
                                Km=coef(nci_out)[2],
                                Ks=NA, 
                                Ki=coef(nci_out)[3], 
                                SSres=sum(resid(nci_out)^2), 
                                SStot=with(data[data$D==i, ], sum((Time-mean(Time))^2)), 
                                P=3, n=nrow(data[data$D==i, ])))
    }, error = function(e){print("NLS cannot converge - trying Differential Optimization algorithm")},
    finally = {
      d = data[data$D==i, ]
      nci_out<-DEoptim(lower=sncil, upper=snciu, fn=ncic,
                       control = c(itermax = 10000, 
                                   trace=FALSE, strategy=3, NP=250))
      sep_out<-rbind(sep_out,
                     data.frame(Model=c("NCI"),
                                Dilution=i, 
                                Vmax=as.numeric(nci_out$optim$bestmem[1]), 
                                Km=as.numeric(nci_out$optim$bestmem[2]),
                                Ks=NA, 
                                Ki=as.numeric(nci_out$optim$bestmem[3]), 
                                SSres=as.numeric(nci_out$optim$bestval), 
                                SStot=with(data[data$D==i, ], sum((Time-mean(Time))^2)), 
                                P=3, n=nrow(data[data$D==i, ])))
    })
    
    ###Competitive inhibition
    tryCatch({
      cia<-nls2(ci, data=data[data$D==i, ], start=snci, control=list(maxiter=1e5), algorithm="brute-force")
      ci_out<-nls2(ci, data=data[data$D==i, ], start=cia, control=list(maxiter=1e5, tol=1e-2, minFactor=1e-12))
      sep_out<-rbind(sep_out,
                     data.frame(Model=c("CI"),
                                Dilution=i, 
                                Vmax=coef(ci_out)[1], 
                                Km=coef(ci_out)[2],
                                Ks=NA, 
                                Ki=coef(ci_out)[3], 
                                SSres=sum(resid(ci_out)^2), 
                                SStot=with(data[data$D==i, ], sum((Time-mean(Time))^2)), 
                                P=3, n=nrow(data[data$D==i, ])))
    }, error = function(e){print("NLS cannot converge - trying Differential Optimization algorithm")},
    finally = {
      d = data[data$D==i, ]
      ci_out<-DEoptim(lower=sncil, upper=snciu, fn=cic,
                       control = c(itermax = 10000, 
                                   trace=FALSE, strategy=3, NP=250))
      sep_out<-rbind(sep_out,
                     data.frame(Model=c("CI"),
                                Dilution=i, 
                                Vmax=as.numeric(ci_out$optim$bestmem[1]), 
                                Km=as.numeric(ci_out$optim$bestmem[2]),
                                Ks=NA, 
                                Ki=as.numeric(ci_out$optim$bestmem[3]), 
                                SSres=as.numeric(ci_out$optim$bestval), 
                                SStot=with(data[data$D==i, ], sum((Time-mean(Time))^2)), 
                                P=3, n=nrow(data[data$D==i, ])))
      
    })
    
    ###Uncompetitive inhibition
    tryCatch({
      ucia<-nls2(uci, data=data[data$D==i, ], start=snci, control=list(maxiter=1e5), algorithm="brute-force")
      uci_out<-nls2(uci, data=data[data$D==i, ], start=ucia, control=list(maxiter=1e5, tol=1e-2, minFactor=1e-12))
      sep_out<-rbind(sep_out,
                     data.frame(Model=c("UCI"),
                                Dilution=i, 
                                Vmax=coef(uci_out)[1], 
                                Km=coef(uci_out)[2],
                                Ks=NA, 
                                Ki=coef(uci_out)[3], 
                                SSres=sum(resid(uci_out)^2), 
                                SStot=with(data[data$D==i, ], sum((Time-mean(Time))^2)), 
                                P=3, n=nrow(data[data$D==i, ])))
    }, error = function(e){print("NLS cannot converge - trying Differential Optimization algorithm")},
    finally = {
      d = data[data$D==i, ]
      uci_out<-DEoptim(lower=sncil, upper=snciu, fn=ucic,
                      control = c(itermax = 10000, 
                                  trace=FALSE, strategy=3, NP=250))
      sep_out<-rbind(sep_out,
                     data.frame(Model=c("UCI"),
                                Dilution=i, 
                                Vmax=as.numeric(uci_out$optim$bestmem[1]), 
                                Km=as.numeric(uci_out$optim$bestmem[2]),
                                Ks=NA, 
                                Ki=as.numeric(uci_out$optim$bestmem[3]), 
                                SSres=as.numeric(uci_out$optim$bestval), 
                                SStot=with(data[data$D==i, ], sum((Time-mean(Time))^2)), 
                                P=3, n=nrow(data[data$D==i, ])))
    })
  }
  
  ##calculate Rsq
  sep_out$Rsq<-1-sep_out$SSres/sep_out$SStot
  
  ##calculate F statitic and p values
  sep_out$Fstat<-numeric(length = nrow(sep_out))
  sep_out$pval<-numeric(length = nrow(sep_out))
  
  ###Dilution 1
  sep_out1<-sep_out[sep_out$Dilution==1, ]
  for(v in 2:nrow(sep_out1)){
    sep_out1[1, "Fstat"]<-NA
    sep_out1[1, "pval"]<-NA
    sep_out1[v, "Fstat"]<-(as.numeric(sep_out1[1, "SSres"])-as.numeric(sep_out1[v, "SSres"]))*
      (as.numeric(sep_out1[1, "n"])-as.numeric(sep_out1[v, "P"]))/
      as.numeric(sep_out1[v, "SSres"])/(as.numeric(sep_out1[v, "P"])-as.numeric(sep_out1[1, "P"]))
    sep_out1[v, "pval"]<-pf(q=((as.numeric(sep_out1[1, "SSres"])-as.numeric(sep_out1[v, "SSres"]))*
                                (as.numeric(sep_out1[1, "n"])-as.numeric(sep_out1[v, "P"]))/
                                as.numeric(sep_out1[v, "SSres"])/(as.numeric(sep_out1[v, "P"])-as.numeric(sep_out1[1, "P"]))), 
                           df1=(as.numeric(sep_out1[v, "P"])-as.numeric(sep_out1[1, "P"])), 
                           df2=(as.numeric(sep_out1[1, "n"])-as.numeric(sep_out1[v, "P"])), 
                           lower.tail=F)
    
    
  }
  ###Dilution 2
  sep_out2<-sep_out[sep_out$Dilution==2, ]
  for(v in 2:nrow(sep_out2)){
    sep_out2[1, "Fstat"]<-NA
    sep_out2[1, "pval"]<-NA
    sep_out2[v, "Fstat"]<-(as.numeric(sep_out2[1, "SSres"])-as.numeric(sep_out2[v, "SSres"]))*
      (as.numeric(sep_out2[1, "n"])-as.numeric(sep_out2[v, "P"]))/
      as.numeric(sep_out2[v, "SSres"])/(as.numeric(sep_out2[v, "P"])-as.numeric(sep_out2[1, "P"]))
    sep_out2[v, "pval"]<-pf(q=((as.numeric(sep_out2[1, "SSres"])-as.numeric(sep_out2[v, "SSres"]))*
                                 (as.numeric(sep_out2[1, "n"])-as.numeric(sep_out2[v, "P"]))/
                                 as.numeric(sep_out2[v, "SSres"])/(as.numeric(sep_out2[v, "P"])-as.numeric(sep_out2[1, "P"]))), 
                            df1=(as.numeric(sep_out2[v, "P"])-as.numeric(sep_out2[1, "P"])), 
                            df2=(as.numeric(sep_out2[1, "n"])-as.numeric(sep_out2[v, "P"])), 
                            lower.tail=F)
    
    
  }
  ###Dilution 4
  sep_out4<-sep_out[sep_out$Dilution==4, ]
  for(v in 2:nrow(sep_out4)){
    sep_out4[1, "Fstat"]<-NA
    sep_out4[1, "pval"]<-NA
    sep_out4[v, "Fstat"]<-(as.numeric(sep_out4[1, "SSres"])-as.numeric(sep_out4[v, "SSres"]))*
      (as.numeric(sep_out4[1, "n"])-as.numeric(sep_out4[v, "P"]))/
      as.numeric(sep_out4[v, "SSres"])/(as.numeric(sep_out4[v, "P"])-as.numeric(sep_out4[1, "P"]))
    sep_out4[v, "pval"]<-pf(q=((as.numeric(sep_out4[1, "SSres"])-as.numeric(sep_out4[v, "SSres"]))*
                                 (as.numeric(sep_out4[1, "n"])-as.numeric(sep_out4[v, "P"]))/
                                 as.numeric(sep_out4[v, "SSres"])/(as.numeric(sep_out4[v, "P"])-as.numeric(sep_out4[1, "P"]))), 
                            df1=(as.numeric(sep_out4[v, "P"])-as.numeric(sep_out4[1, "P"])), 
                            df2=(as.numeric(sep_out4[1, "n"])-as.numeric(sep_out4[v, "P"])), 
                            lower.tail=F)
    
    
  }
  
  sep_outf<-rbind(sep_out1, sep_out2, sep_out4)
  
  return(sep_outf)
  
}