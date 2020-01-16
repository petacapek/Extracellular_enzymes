C_calc<-function(dataset, MUBconc, APconc, Nmeasure, Times, intens){#
  
  #Loop function
  lfunction<-function(rows){
    
    #Reading data
    e_0<-read.csv(file = dataset,
                   header = T, skip = rows[1]-1, nrows = rows[2]-1)
    
    #Rearanging the data
    e_0r<-data.frame(Slurry = c(rep(0, times = 9),# calibration line is listed first
                                rep(3, times = 9),
                                rep(2, times = 9),
                                rep(1, times = 9),
                                rep(0, times = 6),
                                rep(3, times = 3*6),#then, the measurements
                                rep(2, times = 3*6),
                                rep(1, times = 3*6)),
                     Calibration = c(rep("TRUE", times = 12+4*6), #define the data set for calibration and measurment
                                     rep("FALSE", times = 12*5)),
                     C_MUB = c(rep(c(rep(0, times=3), MUBconc/5), times = 4), #concentration of MUB
                               rep(NA, times = 12*5)),
                     C_AP = c(rep(NA, times = 12*3),
                              APconc/5, #concentrations of added AP substrate
                              rep(c(rep(APconc/5, each = 3)), times = 3)),
                     E = c(as.numeric(e_0[1, c(1:3)]), as.numeric(e_0[2, c(1, 5, 9)]),
                           as.numeric(e_0[3, c(1, 5, 9)]),#calibration for extractant
                           as.numeric(e_0[1, c(4:6)]), as.numeric(e_0[2, c(2, 6, 10)]),
                           as.numeric(e_0[3, c(2, 6, 10)]), #calibration for SL 3
                           as.numeric(e_0[1, c(7:9)]), as.numeric(e_0[2, c(3, 7, 11)]),
                           as.numeric(e_0[3, c(3, 7, 11)]), #calibration for SL 2
                           as.numeric(e_0[1, c(10:12)]), as.numeric(e_0[2, c(4, 8, 12)]),
                           as.numeric(e_0[3, c(4, 8, 12)]), #calibration for SL 1
                           as.numeric(e_0[4, c(1:6)]), #data for extractant
                           as.numeric(e_0[4, c(7:9)]), as.numeric(e_0[5, c(4:6)]), 
                           as.numeric(e_0[6, c(1:3)]), as.numeric(e_0[6, c(10:12)]),
                           as.numeric(e_0[7, c(7:9)]), as.numeric(e_0[8, c(4:6)]), #data for SL3 
                           as.numeric(e_0[4, c(10:12)]), as.numeric(e_0[5, c(7:9)]), 
                           as.numeric(e_0[6, c(4:6)]), as.numeric(e_0[7, c(1:3)]),
                           as.numeric(e_0[7, c(10:12)]), as.numeric(e_0[8, c(7:9)]), #data for SL2
                           as.numeric(e_0[5, c(1:3)]), as.numeric(e_0[5, c(10:12)]), 
                           as.numeric(e_0[6, c(7:9)]), as.numeric(e_0[7, c(4:6)]),
                           as.numeric(e_0[8, c(1:3)]), as.numeric(e_0[8, c(10:12)]) #data for SL1
                     ))
    return(e_0r)
  }
  
  #Go through the loop
  ID=1
  for(i in 2:Nmeasure){ID<-append(ID, 1+i*10-10)}
  
  #first measurement
  all<-lfunction(rows = c(1,(1+8)))
  
  #all other measurements
  for(i in ID[-1]){
    all<-rbind(all, lfunction(rows = c(i,(i+8))))
  }
  
  #divide data set into calibrations, blanks and real data
  cals<-subset(all, Calibration == "TRUE")
  blanks<-subset(all, Calibration != "TRUE" & Slurry==0)
  d<-subset(all, Calibration != "TRUE" & Slurry>0)
  
  #Subtract blanks
  cals$Ecorr<-numeric(length = nrow(cals))
  blanks$Ecorr<-numeric(length = nrow(blanks))
  d$Ecorr<-numeric(length = nrow(d))
  ##Calibrations
  for(i in unique(cals$Slurry)){
    cals[cals$Slurry==i, "Ecorr"]<-cals[cals$Slurry==i, "E"]-
      mean(cals[(cals$Slurry==i & cals$C_MUB==0), "E"], na.rm=T)
  }
  ##blanks
  blanks$Ecorr<-blanks$E-mean(cals[(cals$Slurry==0 & cals$C_MUB==0), "E"], na.rm=T)
  ##data
  d$Time<-Times
  
  if(intens=="TRUE"){
    for(i in unique(d$Slurry)){
      d[d$Slurry==i, "Ecorr"]<-d[d$Slurry==i, "E"]-
        mean(c(as.numeric(cals[(cals$Slurry==i & cals$C_MUB==0), "E"]),
               as.numeric(d[(d$Slurry==i & d$Time==0), "E"])), na.rm=T)
    }
  }else{
    for(i in unique(d$Slurry)){
      d[d$Slurry==i, "Ecorr"]<-d[d$Slurry==i, "E"]-
        mean(as.numeric(cals[(cals$Slurry==i & cals$C_MUB==0), "E"]), na.rm=T)
    }
  }
  
  
  #Use linear regression to make a calibration line for each dilution separately
  cal0_0<-lm(Ecorr~C_MUB-1, data = cals[(cals$Slurry==0), ])
  cal0_1<-lm(Ecorr~C_MUB-1, data = cals[(cals$Slurry==1), ])
  cal0_2<-lm(Ecorr~C_MUB-1, data = cals[(cals$Slurry==2), ])
  cal0_3<-lm(Ecorr~C_MUB-1, data = cals[(cals$Slurry==3), ])
  #extract coefficients
  cal0_0_coef<-coef(cal0_0)
  cal0_1_coef<-coef(cal0_1)
  cal0_2_coef<-coef(cal0_2)
  cal0_3_coef<-coef(cal0_3)
  
  #calculate concentration of product according to regressions
  #concentration of product is calculated in umols of product per l
  ##data
  d$P<-numeric(length = nrow(d))
  for(i in 1:nrow(d)){
    if(d$Slurry[i]==1){
      d$P[i]<-d$Ecorr[i]/cal0_1_coef
      }else{
        if(d$Slurry[i]==2){
          d$P[i]<-d$Ecorr[i]/cal0_2_coef
          }else{
            d$P[i]<-d$Ecorr[i]/cal0_3_coef
          }
        }
      }
  ##blanks
  blanks$P<-blanks$Ecorr/cal0_0_coef
  
  #Add blanks to data
  d$Pblanks<-numeric(length = nrow(d))
  for(i in unique(d$C_AP)){
    d[d$C_AP==i, "Pblanks"]<-rep(mean(blanks[blanks$C_AP==i, "P"], na.rm=T), 
                                 times=nrow(d[d$C_AP==i, ]))
  }
  
  #Subtract the blanks
  d$Pcorr<-d$P-d$Pblanks
  
  
  all_out<-list(data=d, 
                blanks=blanks, 
                Cals=list(S0=cal0_0, S1=cal0_1, S2=cal0_2, S3=cal0_3), 
                Coefs=c(S0=cal0_0_coef, S1=cal0_1_coef, S2=cal0_2_coef, S3=cal0_3_coef),
                PlotS0=(ggplot(data=cals[(cals$Slurry==0), ], aes(C_MUB, Ecorr))+
                  geom_point(cex=6, pch=21, fill="grey")+
                  stat_smooth(method = lm, se=F)),
                PlotS1=(ggplot(data=cals[(cals$Slurry==1), ], aes(C_MUB, Ecorr))+
                  geom_point(cex=6, pch=21, fill="grey")+
                  stat_smooth(method = lm, se=F)),
                PlotS2=(ggplot(data=cals[(cals$Slurry==2), ], aes(C_MUB, Ecorr))+
                  geom_point(cex=6, pch=21, fill="grey")+
                  stat_smooth(method = lm, se=F)),
                PlotS3=(ggplot(data=cals[(cals$Slurry==3), ], aes(C_MUB, Ecorr))+
                  geom_point(cex=6, pch=21, fill="grey")+
                  stat_smooth(method = lm, se=F))
                )
  return(all_out)
  
}