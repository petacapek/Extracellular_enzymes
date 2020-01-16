#LIBRARIES
library(ggplot2)
library(deSolve)
library(FME)
library(DEoptim)
library(minpack.lm)
library(soilphysics)
library(dplyr)
library(lamW)
library(nls2)
library(segmented)
library(reshape)
library(gridExtra)

##GGPLOT THEME
theme_min<-theme(axis.text.x=element_text(vjust=0.2, size=18, colour="black"),
                 axis.text.y=element_text(hjust=0.2, size=18, colour="black"),
                 axis.title=element_text(size=18, colour="black"),
                 axis.line=element_line(size=0.5, colour="black"),
                 strip.text=element_text(size=18, face="bold"),
                 axis.ticks=element_line(size=1, colour="black"),
                 axis.ticks.length=unit(-0.05, "cm"),
                 panel.background=element_rect(colour="black", fill="white"),
                 panel.grid=element_line(linetype=0),
                 legend.text=element_text(size=14, colour="black"),
                 legend.title=element_text(size=14, colour="black"),
                 legend.position=c("right"),
                 legend.key.size=unit(1, "cm"),
                 strip.background=element_rect(fill="grey98", colour="black"),
                 legend.key=element_rect(fill="white", size=1.2),
                 legend.spacing=unit(0.5, "cm"),
                 plot.title=element_text(size=18, face="bold", hjust=-0.05))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#####################################Product concentration################################
source("./r_functions/C_calc.R")
##########################################Assay 1##################################################
#Calculating data
exp1_data<-C_calc(dataset = c("./data_files/5.09.19_assay_no1/assay1_all.csv"),
           MUBconc = c(1, 5, 10, 25, 50, 100),
           APconc = c(10, 25, 50, 100, 250, 500),
           Nmeasure = 6,
           Times = rep(c(0, 40/60, 70/60, 127/60, 247/60, 365/60),
                       each = 3*3*6),
           intens = "FALSE")
#Export product formation data
e1<-exp1_data$data
e1$Pcorr2<-e1$Pcorr
#Set negative numbers to zero
e1[e1$Pcorr2<0, "Pcorr2"]<-0
#Define soil - buffer solution dilution rate
e1$D<-e1$Slurry
e1[e1$D==3, "D"]<-4

#Check the data 
ggplot(e1, aes(Time, Pcorr2))+geom_point(cex=6, pch=21, aes(fill=factor(D)))+
  facet_wrap(D~C_AP, scales="free")
#Mark the outliers
e1$Outlier<-c(0)
e1[(e1$C_AP==5 & e1$D==1 & e1$Time>1 & e1$Time<2 & e1$Pcorr2<0.5), "Outlier"]<-c(1)
e1[(e1$C_AP==5 & e1$D==4 & e1$Time>1 & e1$Time<2 & e1$Pcorr2>0.5), "Outlier"]<-c(1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##########################################Assay 2##################################################
#Calculating data
exp2_data<-C_calc(dataset = c("./data_files/5.10.19_assay_no2/assay2_all.csv"),
                  Times = rep(c(0, 0.5, 1, 117/60, 4, 6), each=3*6*3),
                  MUBconc = c(1, 5, 10, 25, 50, 100),
                  APconc = c(10, 25, 50, 100, 250, 500),
                  Nmeasure = 6,
                  intens = "FALSE")
#Export product formation data
e2<-exp2_data$data
e2$Pcorr2<-e2$Pcorr
#Set negative numbers to zero
e2[e2$Pcorr2<0, "Pcorr2"]<-0
#Define soil - buffer solution dilution rate
e2$D<-e2$Slurry
e2[e2$D==3, "D"]<-4

#Check the data 
ggplot(e2, aes(Time, Pcorr2))+geom_point(cex=6, pch=21, aes(fill=factor(D)))+
  facet_wrap(D~C_AP, scales="free")
#Mark the outliers
e2$Outlier<-c(0)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##########################################Assay 3##################################################
#Calculating data
exp3_data<-C_calc(dataset = c("./data_files/5.16.19_assay_no3/assay3_all.csv"),
                  MUBconc = c(0.77, 3.85, 7.7, 23.1, 46.2, 92.4),
                  APconc = c(6.98, 20.94, 41.88, 83.75, 251.26, 502.52),
                  Nmeasure = 6,
                  Times = rep(c(0, 0.5, 1, 2, 4, 6), each=3*6*3), intens = "FALSE")
#Export product formation data
e3<-exp3_data$data
e3$Pcorr2<-e3$Pcorr
#Set negative numbers to zero
e3[e3$Pcorr2<0, "Pcorr2"]<-0
#Define soil - buffer solution dilution rate
e3$D<-e3$Slurry
e3[e3$D==3, "D"]<-4

#Check the data 
ggplot(e3, aes(Time, Pcorr2))+geom_point(cex=6, pch=21, aes(fill=factor(D)))+
  facet_wrap(D~C_AP, scales="free")
#Mark the outliers
e3$Outlier<-c(0)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##########################################Assay 4##################################################
#Calculating data
exp4_data<-C_calc(dataset = c("./data_files/5.24.19_assay_no4/assay4_all.csv"),
                  Times = rep(c(0, 0.5, 1, (120-17)/60, 4, 6), each=3*6*3),
                  MUBconc = c(0.757, 3.784, 7.568, 22.705, 45.41, 90.82),
                  APconc = c(10.30, 51.51, 103.018, 257.544, 515.089, 1030.178),
                  Nmeasure = 6,
                  intens = "FALSE")
#Export product formation data
e4<-exp4_data$data
e4$Pcorr2<-e4$Pcorr
#Set negative numbers to zero
e4[e4$Pcorr2<0, "Pcorr2"]<-0
#Define soil - buffer solution dilution rate
e4$D<-e4$Slurry
e4[e4$D==3, "D"]<-4

#Check the data 
ggplot(e4, aes(Time, Pcorr2))+geom_point(cex=6, pch=21, aes(fill=factor(D)))+
  facet_wrap(D~C_AP, scales="free")
#Mark the outliers
e4$Outlier<-c(0)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#####################################Assay 5####################################################
#Calculating data
#Time
t5=read.csv(file = c("./data_files/5.30.19_assay_no5/assay5_t.csv"))
a5_t<-numeric()
for(i in 1:19){
  a5_t<-append(a5_t, c(rep(rep(as.numeric(t5[t5$Measurement==i, "Time"]), each=3),
                               times = 3)))
}

exp5_data<-C_calc(dataset = c("./data_files/5.30.19_assay_no5/assay5_all.csv"),
           MUBconc = c(0.7956, 3.9782, 7.95633, 23.86899, 47.73798, 95.47596),
           APconc = c(10.356, 51.78216, 103.5643, 258.9108, 517.8216, 1035.643),
           Nmeasure = 19,
           Times = a5_t, intens = "TRUE")
#Export product formation data
e5<-exp5_data$data
e5$Pcorr2<-e5$Pcorr
#Set negative numbers to zero
e5[e5$Pcorr2<0, "Pcorr2"]<-0
#Define soil - buffer solution dilution rate
e5$D<-e5$Slurry
e5[e5$D==3, "D"]<-4

#Check the data 
ggplot(e5, aes(Time, Pcorr2))+geom_point(cex=6, pch=21, aes(fill=factor(D)))+
  facet_wrap(D~C_AP, scales="free")
#Mark the outliers
e5$Outlier<-c(0)

e5[(e5$D==2 & e5$C_AP==unique(e5$C_AP)[4] & e5$Pcorr2>3), "Outlier"]<-c(1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##########################################Assay 6##################################################
#Calculating data
#Times
t6<-read.csv(file=c("./data_files/7.10.19_assay_no6/assay6_t.csv"))
t6<-t6[order(t6$Measurement, -t6$Slurry), ]
a6_t<-numeric()

for(i in 1:max(t6$Measurement)){
  a6_t<-append(a6_t, rep(as.numeric(t6[t6$Measurement==i, "Time"]), each=3))
}

#Calculating data
exp6_data<-C_calc(dataset = "./data_files/7.10.19_assay_no6/assay6_all.csv",
                  Times = a6_t,
                  MUBconc = c(1.83, 4.575, 9.15, 22.8756, 45.751, 91.503),
                  APconc = c(13.539, 33.8474, 67.6947, 135.3894, 228.474, 1015.42),
                  Nmeasure = max(t6$Measurement),
                  intens = "TRUE")
#Export product formation data
e6<-exp6_data$data
e6$Pcorr2<-e6$Pcorr
#Set negative numbers to zero
e6[e6$Pcorr2<0, "Pcorr2"]<-0
#Define soil - buffer solution dilution rate
e6$D<-e6$Slurry
e6[e6$D==3, "D"]<-4

#Check the data 
ggplot(e6, aes(Time, Pcorr2))+geom_point(cex=6, pch=21, aes(fill=factor(D)))+
  facet_wrap(D~C_AP, scales="free")
#Mark the outliers
e6$Outlier<-c(0)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##############################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###############################Non-linear analysis - first 4 assays###########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
source("./r_functions/Emultinls.R")

#Assay 1
e1_nls<-Emultinls(data=e1[e1$Outlier==0, ])
##export the results to csv
write.csv(e1_nls, "./e1_nls.csv")

#Assay 2
e2_nls<-Emultinls(data=e2[e2$Outlier==0, ])
##export the results to csv
write.csv(e2_nls, "./e2_nls.csv")

#Assay 3
e3_nls<-Emultinls(data=e3[e3$Outlier==0, ])
##export the results to csv
write.csv(e3_nls, "./e3_nls.csv")

#Assay 4
e4_nls<-Emultinls(data=e4[e4$Outlier==0, ])
##export the results to csv
write.csv(e4_nls, "./e4_nls.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###############################Non-linear analysis - last 2 assays############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Define models
##Two enzymes - Non-competitive inhibition of one
twoe_NCI<-function(time, state, pars){
  
  with(as.list(c(state, pars)),{
    
    dP<-(Vmax1*S/(Km + S))+(Vmax2*S/(Km*(1+P/Kic) + S*(1+P/Kiu)))
    dS<--(Vmax1*S/(Km + S))-(Vmax2*S/(Km*(1+P/Kic) + S*(1+P/Kiu)))
    
    return(list(c(dP, dS)))
    
  })
}
##Two enzymes - Competitive inhibition of one
twoe_CI<-function(time, state, pars){
  
  with(as.list(c(state, pars)),{
    
    dP<-(Vmax1*S/(Km + S))+(Vmax2*S/(Km*(1+P/Kic) + S))
    dS<--(Vmax1*S/(Km + S))-(Vmax2*S/(Km*(1+P/Kic) + S))
    
    return(list(c(dP, dS)))
    
  })
}
##Two enzymes - Un-competitive inhibition of one
twoe_UCI<-function(time, state, pars){
  
  with(as.list(c(state, pars)),{
    
    dP<-(Vmax1*S/(Km + S))+(Vmax2*S/(Km + S*(1+P/Kic)))
    dS<--(Vmax1*S/(Km + S))-(Vmax2*S/(Km + S*(1+P/Kic)))
    
    return(list(c(dP, dS)))
    
  })
}

#Load the estimation functions
source("./r_functions/twoeOpt_NCI.R")
source("./r_functions/twoeOpt.R")

#Assay 5
#Single enzyme - substrate pair - non-inhibited Michaelis-Menten kinetic
e5d1MM<-nlsLM(Pcorr2~C_AP-Pblanks-Km*lambertW0((C_AP-Pblanks)/Km*exp((C_AP-Pblanks-Vmax*Time)/Km)), 
              data = subset(e5, D==1 & Outlier==0),
              start=list(Vmax=1, Km=10))
summary(e5d1MM)
Rsq(e5d1MM)
sum(resid(e5d1MM)^2)

e5d2MM<-nlsLM(Pcorr2~C_AP-Pblanks-Km*lambertW0((C_AP-Pblanks)/Km*exp((C_AP-Pblanks-Vmax*Time)/Km)), 
              data = subset(e5, D==2 & Outlier==0),
              start=list(Vmax=1, Km=10))
summary(e5d2MM)
Rsq(e5d2MM)
sum(resid(e5d2MM)^2)

e5d4MM<-nlsLM(Pcorr2~C_AP-Pblanks-Km*lambertW0((C_AP-Pblanks)/Km*exp((C_AP-Pblanks-Vmax*Time)/Km)), 
              data = subset(e5, D==4 & Outlier==0),
              start=list(Vmax=1, Km=10))
summary(e5d4MM)
Rsq(e5d4MM)
sum(resid(e5d4MM)^2)

#CI
e5d1_twoe<-twoeOpt(data=e5[(e5$D==1 & e5$Outlier==0), ], model=twoe_CI)
e5d1_twoe$P
e5d1_twoe$Goodness$Gfit

e5d2_twoe<-twoeOpt(data=e5[(e5$D==2 & e5$Outlier==0), ], model=twoe_CI)
e5d2_twoe$P
e5d2_twoe$Goodness$Gfit

e5d4_twoe<-twoeOpt(data=e5[(e5$D==4 & e5$Outlier==0), ], model=twoe_CI)
e5d4_twoe$P
e5d4_twoe$Goodness$Gfit

#UCI
e5d1_twoeb<-twoeOpt(data=e5[(e5$D==1 & e5$Outlier==0), ], model=twoe_UCI)
e5d1_twoeb$P
e5d1_twoeb$Goodness$Gfit

e5d2_twoeb<-twoeOpt(data=e5[(e5$D==2 & e5$Outlier==0), ], model=twoe_UCI)
e5d2_twoeb$P
e5d2_twoeb$Goodness$Gfit

e5d4_twoeb<-twoeOpt(data=e5[(e5$D==4 & e5$Outlier==0), ], model=twoe_UCI)
e5d4_twoeb$P
e5d4_twoeb$Goodness$Gfit

#NCI
e5d1_twoec<-twoeOpt_NCI(data=e5[(e5$D==1 & e5$Outlier==0), ], model=twoe_NCI)
e5d1_twoec$P
e5d1_twoec$Goodness$Gfit

e5d2_twoec<-twoeOpt_NCI(data=e5[(e5$D==2 & e5$Outlier==0), ], model=twoe_NCI)
e5d2_twoec$P
e5d2_twoec$Goodness$Gfit

e5d4_twoec<-twoeOpt_NCI(data=e5[(e5$D==4 & e5$Outlier==0), ], model=twoe_NCI)
e5d4_twoec$P
e5d4_twoec$Goodness$Gfit

#Assay 6
#Single enzyme - substrate pair - non-inhibited Michaelis-Menten kinetic
e6d1MM<-nlsLM(Pcorr2~C_AP-Pblanks-Km*lambertW0((C_AP-Pblanks)/Km*exp((C_AP-Pblanks-Vmax*Time)/Km)), 
              data = subset(e6, D==1 & Outlier==0),
              start=list(Vmax=1, Km=10))
summary(e6d1MM)
Rsq(e6d1MM)
sum(resid(e6d1MM)^2)

e6d2MM<-nlsLM(Pcorr2~C_AP-Pblanks-Km*lambertW0((C_AP-Pblanks)/Km*exp((C_AP-Pblanks-Vmax*Time)/Km)), 
              data = subset(e6, D==2 & Outlier==0),
              start=list(Vmax=1, Km=10))
summary(e6d2MM)
Rsq(e6d2MM)
sum(resid(e6d2MM)^2)

e6d4MM<-nlsLM(Pcorr2~C_AP-Pblanks-Km*lambertW0((C_AP-Pblanks)/Km*exp((C_AP-Pblanks-Vmax*Time)/Km)), 
              data = subset(e6, D==4 & Outlier==0),
              start=list(Vmax=1, Km=10))
summary(e6d4MM)
Rsq(e6d4MM)
sum(resid(e6d4MM)^2)

#CI
e6d1_twoe<-twoeOpt(data=e6[(e6$D==1 & e6$Outlier==0), ], model=twoe_CI)
e6d1_twoe$P
e6d1_twoe$Goodness$Gfit

e6d2_twoe<-twoeOpt(data=e6[(e6$D==2 & e6$Outlier==0), ], model=twoe_CI)
e6d2_twoe$P
e6d2_twoe$Goodness$Gfit

e6d4_twoe<-twoeOpt(data=e6[(e6$D==4 & e6$Outlier==0), ], model=twoe_CI)
e6d4_twoe$P
e6d4_twoe$Goodness$Gfit

#UCI
e6d1_twoeb<-twoeOpt(data=e6[(e6$D==1 & e6$Outlier==0), ], model=twoe_UCI)
e6d1_twoeb$P
e6d1_twoeb$Goodness$Gfit

e6d2_twoeb<-twoeOpt(data=e6[(e6$D==2 & e6$Outlier==0), ], model=twoe_UCI)
e6d2_twoeb$P
e6d2_twoeb$Goodness$Gfit

e6d4_twoeb<-twoeOpt(data=e6[(e6$D==4 & e6$Outlier==0), ], model=twoe_UCI)
e6d4_twoeb$P
e6d4_twoeb$Goodness$Gfit

#NCI
e6d1_twoec<-twoeOpt_NCI(data=e6[(e6$D==1 & e6$Outlier==0), ], model=twoe_NCI)
e6d1_twoec$P
e6d1_twoec$Goodness$Gfit

e6d2_twoec<-twoeOpt_NCI(data=e6[(e6$D==2 & e6$Outlier==0), ], model=twoe_NCI)
e6d2_twoec$P
e6d2_twoec$Goodness$Gfit

e6d4_twoec<-twoeOpt_NCI(data=e6[(e6$D==4 & e6$Outlier==0), ], model=twoe_NCI)
e6d4_twoec$P
e6d4_twoec$Goodness$Gfit

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################Figures##########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Main text
#Fig. 1
Theor2<-data.frame(t=seq(0, 6, by=0.1))
Theor2$P<-10-5*lambertW0(10/5*exp((10-10*Theor2$t)/5))
Theor2$Legend<-c("10")
Theor2<-rbind(Theor2, data.frame(t=seq(0, 6, by=0.1),
                               P=10-5*lambertW0(10/5*exp((10-5*Theor2$t)/5)),
                               Legend=c("5")))
Theor2<-rbind(Theor2, data.frame(t=seq(0, 6, by=0.1),
                               P=10-5*lambertW0(10/5*exp((10-0.5*Theor2$t)/5)),
                               Legend=c("1")))

Theor2$Legend<-factor(Theor2$Legend, levels = c("10", "5", "1"))

(TA<-ggplot(Theor2, aes(t, P))+geom_line(lwd=1, aes(linetype=Legend), colour="black")+
  ylim(0, 10)+theme_min+theme(legend.title = element_blank(),
                              legend.position = c(0.7, 0.4),
                              legend.text.align = 0)+
  xlab("Time (h)")+scale_linetype_manual(values = c("solid", "dashed", "dotdash"),
                                         labels = expression(V[MAX]==10~mu*mol~L^{-1}~h^{-1},
                                                             V[MAX]==5~mu*mol~L^{-1}~h^{-1},
                                                             V[MAX]==1~mu*mol~L^{-1}~h^{-1}))+
  ylab(expression(paste("Product (", mu, "mol ", L^{-1}, ")")))+
  ggtitle("A)"))


##Inhibition
###Without inhibition
WIout<-data.frame(time=seq(0, 6, by=0.1),
           P=10-5*lambertW0(10/5*exp((10-1*Theor2$t)/5)),
           Legend=c("WI"))

###Substrate inhibition
SI<-function(time, state, pars){
  
  with(as.list(c(state, pars)),{
    
    dP<-Vmax*(10-P)/(Km + (10-P) + (10-P)^2/Ki)
    
    return(list(c(dP)))
    
  })
}

SIout <- as.data.frame(ode(y=c(P=0), parms=c(Vmax=1, Km=5, Ki=5), 
             SI, times=seq(0,6, by=0.1)))
plot(SIout)

SIout$Legend<-"SI"

ggplot(SIout, aes(time, P)) + geom_point(cex=2) + stat_smooth(method=lm)

###Noncompetitive
NCI<-function(time, state, pars){
  
  with(as.list(c(state, pars)),{
    
    dP<-Vmax*(10-P)/(Km*(1+P/Kic) + (10-P)*(1+P/Kiu))
    
    return(list(c(dP)))
    
  })
}

NCIout <- as.data.frame(ode(y=c(P=0), parms=c(Vmax=1, Km=5, Kic=1, Kiu=1), 
             NCI, times=seq(0,6, by=0.1)))
plot(NCIout)

NCIout$Legend<-"NCI"

###Competitive
CI<-function(time, state, pars){
  
  with(as.list(c(state, pars)),{
    
    dP<-Vmax*(10-P)/(Km*(1+P/Kic) + (10-P))
    
    return(list(c(dP)))
    
  })
}

CIout <- as.data.frame(ode(y=c(P=0), parms=c(Vmax=1, Km=5, Kic=1), 
              CI, times=seq(0,6, by=0.1)))
plot(CIout)

CIout$Legend<-"CI"

###Uncompetitive
UCI<-function(time, state, pars){
  
  with(as.list(c(state, pars)),{
    
    dP<-Vmax*(10-P)/(Km + (10-P)*(1+P/Kiu))
    
    return(list(c(dP)))
    
  })
}

UCIout <- as.data.frame(ode(y=c(P=0), parms=c(Vmax=1, Km=5, Kiu=1), 
              UCI, times=seq(0,6, by=0.1)))
plot(UCIout)

UCIout$Legend<-"UCI"

Theor3<-rbind(WIout, SIout, NCIout, CIout, UCIout)

(TB<-ggplot(Theor3, aes(time, P)) + geom_line(lwd=1, aes(linetype=Legend, colour=Legend)) +
  theme_min + theme(legend.title = element_blank(),
                   legend.position = c(0.2, 0.7)) +
  xlab("Time (h)") + ylab(expression(paste("Product (", mu, "mol ", L^{-1}, ")")))+
  ggtitle("B)") + 
  scale_linetype_manual(values = c("solid", "solid", "dashed", "dashed", 
                                   "dotdash"))+
  scale_color_manual(values = c("black", "grey50", "grey50","black", "black")))

grid.arrange(TA, TB, nrow=1)



#Figure 1: Progress curves of enzymatic reaction following the non-inhibited Michaelis-Menten kinetic. 
#The effect of maximum velocity constant (VMAX) and enzyme affinity to substrate (KM) is depicted by 
#different line types and line colors respectively. For all curves, initial substrate concentration (S0)
#is 10 µmol L-1. Black lines denote progress curves with KM equal to 5 µmol L-1, and solid grey lines 
#curves with KM equal to 10 µmol L-1. 

#Fig. 2
e4$C_AP2<-round(e4$C_AP, 0)
model4<-nls(Pcorr2~(C_AP-Pblanks)-Km[factor(D)]*lambertW0((C_AP-Pblanks)/Km[factor(D)]*exp((C_AP-Pblanks-Vmax[factor(D)]*Time)/Km[factor(D)])), 
            data = e4[e4$Outlier==0, ],
            start=list(Vmax=c(1, 1, 1), Km=c(10, 10, 10)))
summary(model4)
e4_sim<-data.frame(Time = rep(seq(0, 6, by=0.1), times=6),
                   C_AP = rep(rep(unique(e4$C_AP), each=length(seq(0, 6, by=0.1))), times = 3),
                   Pblanks = rep(rep(unique(e4$Pblanks), each=length(seq(0, 6, by=0.1))), times = 3),
                   D = (rep(unique(e4$D), each = length(seq(0, 6, by=0.1))*6)))
e4_sim$C_AP2<-round(e4_sim$C_AP, 0)
e4_sim$Pred<-predict(model4, newdata = e4_sim)

e4 %>% filter(Outlier==0) %>% group_by(Time, C_AP2, D) %>% 
  summarize(y=mean(Pcorr2, na.rm=T),
            y.sd=sd(Pcorr2, na.rm=T)) %>%
  ggplot(aes(Time, y))+geom_point(cex=6, pch=21, aes(fill=factor(D)))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), lwd=1)+
  scale_fill_manual(values=c("white", "grey", "black"), 
                    labels = c("1:100", "1:200", "1:400"))+
  facet_wrap(~C_AP2, scales="free")+theme_min+theme(legend.title = element_blank())+
  xlab("Time (h)")+ylab(expression(paste("Product (", mu, "mol ", L^{-1}, ")")))+
  geom_line(data=e4_sim, aes(x=Time, y=Pred, color=factor(D)), lwd=1)+
  scale_colour_manual(values = c("grey80", "grey40", "black"), 
                      labels = c("1:100", "1:200", "1:400"))+
  theme(legend.title = element_blank(),
        legend.position = c(0.05, 0.9))

#Figure 2: Progress curves of MUB-P substrate degradation measured on day 14. 
#The figure is divided into six boxes representing six different initial concentrations (S0) 
#of MUB-P substrate (denoted by the box headings). Empty, grey and black symbols denote three different fresh soil - buffer 
#solution slurries. Solid lines represent the fit of integrated form of non-inhibited 
#Michaelis-Menten equation of single enzyme - substrate pair (eq. 4), which corresponded best 
#to measured data according to non-linear least square regression analysis using five different 
#equations reported in Table 1 (see Tab. 2 and section 2.2). Different line colors denote three 
#different fresh soil - buffer solution slurries, whose data has been fitted by non-linear 
#regression separately. Symbols represent means and error bars their standard deviations. 
#Note that boxes have different y-axis scales.

#Fig. 3
e4$res<-resid(model4)
e4 %>% filter(Outlier==0) %>% group_by(Time, C_AP2, D) %>% 
  summarize(y=mean(res, na.rm=T),
            y.sd=sd(res, na.rm=T)) %>%
  ggplot(aes(Time, y))+geom_point(cex=6, pch=21, aes(fill=factor(D)))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), lwd=1)+
  scale_fill_manual(values=c("white", "grey", "black"), 
                    labels = c("1:100", "1:200", "1:400"))+
  facet_wrap(~C_AP2, scales="free")+theme_min+
  theme(legend.title = element_blank(),
        legend.position = c(0.05, 0.9))+
  geom_hline(yintercept = 0)+
  xlab("Time (h)")+ylab(expression(paste("Residuals (", mu, "mol ", L^{-1}, ")")))

#Figure 3: Residuals of the fit of integrated form of non-inhibited Michaelis-Menten equation 
#depicted in Fig. 2. The figure is divided into six boxes representing six different initial 
#concentrations (S0) of MUB-P substrate (denoted by the box headings). Empty, grey and black symbols denote three different 
#fresh soil - buffer solution slurries. Solid horizonal line highlight the value 
#zero - i.e. absolute correspondence between the data and the equation fit. Symbols represent 
#means and error bars their standard deviations. Note that boxes have different y-axis scales.

#Fig. 4
#Linear analysis
##Linear regression without the intercept
e4$vu<-numeric(length = nrow(e4))
e4$v<-numeric(length = nrow(e4))
e4$vu.se<-numeric(length = nrow(e4))
e4$v.se<-numeric(length = nrow(e4))
for(i in unique(e4$C_AP)){
  for(n in 1:3){
    #Run linear regression
    lm_u<-lm(Pcorr2~Time-1, data=subset(e4, C_AP==i & Slurry==n))
    #Extract coeficients
    e4[(e4$C_AP==i & e4$Slurry==n), "vu"]<-rep(coef(lm_u)[1], 
                                               times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
    e4[(e4$C_AP==i & e4$Slurry==n), "vu.se"]<-rep(summary(lm_u)$coefficients[2], 
                                                  times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
    #Extract residuals and times
    res<-resid(lm_u)
    tres<-as.numeric(e4[(e4$C_AP==i & e4$Slurry==n), "Time"])
    lm0<-lm(res~tres)
    #When breakpoint is found
    if(as.numeric(davies.test(lm0)[5])<0.01){
      lm_c<-lm(Pcorr2~Time-1, data=subset(e4, C_AP==i & Slurry==n &
                                            Time<=as.numeric(davies.test(lm0)[3])))
      #Extract coeficients
      e4[(e4$C_AP==i & e4$Slurry==n), "v"]<-rep(coef(lm_c)[1], 
                                                times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
      e4[(e4$C_AP==i & e4$Slurry==n), "v.se"]<-rep(summary(lm_c)$coefficients[2], 
                                                   times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
      
    }else{
      e4[(e4$C_AP==i & e4$Slurry==n), "v"]<-rep(coef(lm_u)[1], 
                                                times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
      e4[(e4$C_AP==i & e4$Slurry==n), "v.se"]<-rep(summary(lm_u)$coefficients[2], 
                                                   times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
    }
  }
}
##Linear regression with the intercept
e4$vuI<-numeric(length = nrow(e4))
e4$vI<-numeric(length = nrow(e4))
e4$vuI.se<-numeric(length = nrow(e4))
e4$vI.se<-numeric(length = nrow(e4))
e4$Ivu<-numeric(length = nrow(e4))
e4$Iv<-numeric(length = nrow(e4))
e4$Ivu.se<-numeric(length = nrow(e4))
e4$Iv.se<-numeric(length = nrow(e4))

for(i in unique(e4$C_AP)){
  for(n in 1:3){
    #Run linear regression
    lm_u<-lm(Pcorr2~Time, data=subset(e4, C_AP==i & Slurry==n))
    #Extract coeficients
    e4[(e4$C_AP==i & e4$Slurry==n), "Ivu"]<-rep(coef(lm_u)[1], 
                                                times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
    e4[(e4$C_AP==i & e4$Slurry==n), "Ivu.se"]<-rep(summary(lm_u)$coefficients[1, 2], 
                                                   times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
    e4[(e4$C_AP==i & e4$Slurry==n), "vuI"]<-rep(coef(lm_u)[2], 
                                                times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
    e4[(e4$C_AP==i & e4$Slurry==n), "vuI.se"]<-rep(summary(lm_u)$coefficients[2, 2], 
                                                   times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
    #Extract residuals and times
    res<-resid(lm_u)
    tres<-as.numeric(e4[(e4$C_AP==i & e4$Slurry==n), "Time"])
    lm0<-lm(res~tres)
    #When breakpoint is found
    if(as.numeric(davies.test(lm0)[5])<0.01){
      lm_c<-lm(Pcorr2~Time, data=subset(e4, C_AP==i & Slurry==n &
                                          Time<=as.numeric(davies.test(lm0)[3])))
      #Extract coeficients
      e4[(e4$C_AP==i & e4$Slurry==n), "Iv"]<-rep(coef(lm_c)[1], 
                                                 times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
      e4[(e4$C_AP==i & e4$Slurry==n), "Iv.se"]<-rep(summary(lm_c)$coefficients[1, 2], 
                                                    times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
      e4[(e4$C_AP==i & e4$Slurry==n), "vI"]<-rep(coef(lm_c)[2], 
                                                 times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
      e4[(e4$C_AP==i & e4$Slurry==n), "vI.se"]<-rep(summary(lm_c)$coefficients[2, 2], 
                                                    times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
    }else{
      e4[(e4$C_AP==i & e4$Slurry==n), "Iv"]<-rep(coef(lm_u)[1], 
                                                 times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
      e4[(e4$C_AP==i & e4$Slurry==n), "Iv.se"]<-rep(summary(lm_u)$coefficients[1, 2], 
                                                    times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
      e4[(e4$C_AP==i & e4$Slurry==n), "vI"]<-rep(coef(lm_u)[2], 
                                                 times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
      e4[(e4$C_AP==i & e4$Slurry==n), "vI.se"]<-rep(summary(lm_u)$coefficients[2, 2], 
                                                    times=nrow(e4[(e4$C_AP==i & e4$Slurry==n), ]))
    }
  }
}

Velocities<-melt(e4[, c("C_AP", "D", "v", "vu", "vI", "vuI")],
                 id.vars = c("C_AP", "D"))
Velocities$se<-melt(e4[, c("C_AP", "D", "v.se", "vu.se", "vI.se", "vuI.se")],
                    id.vars = c("C_AP", "D"))[, 4]

Velocities$fac<-character(length = nrow(Velocities))
for(i in 1:nrow(Velocities)){
  if(Velocities$D[i]==1){
    Velocities$fac[i]<-c("1:100")
  }else{
    if(Velocities$D[i]==2){
      Velocities$fac[i]<-c("1:200")
    }else{
      Velocities$fac[i]<-c("1:400")
    }
  }
}

nlsMM<-nls(Pcorr2~C_AP-Km[factor(D)]*lambertW0(C_AP/Km[factor(D)]*exp((C_AP-Vmax[factor(D)]*Time)/Km[factor(D)])),
           data = e4, start = list(Vmax=c(1,1,1), Km=c(10, 10, 10)))
summary(nls(vu~Vmax[factor(D)]*C_AP/(Km[factor(D)] + C_AP),
            start = list(Vmax=c(1,1,1), Km=c(10, 10, 10)), data=e4))
summary(nls(v~Vmax[factor(D)]*C_AP/(Km[factor(D)] + C_AP),
            start = list(Vmax=c(1,1,1), Km=c(10, 10, 10)), data=e4))
summary(nls(vuI~Vmax[factor(D)]*C_AP/(Km[factor(D)] + C_AP),
            start = list(Vmax=c(1,1,1), Km=c(10, 10, 10)), data=e4))
summary(nls(vI~Vmax[factor(D)]*C_AP/(Km[factor(D)] + C_AP),
            start = list(Vmax=c(1,1,1), Km=c(10, 10, 10)), data=e4))

velsim<-data.frame(C_AP = rep(seq(0, 206, by=1), times=3),
                   fac = rep(c("1:100", "1:200", "1:400"), 
                             each = length(seq(0, 206, by=1))))
velsim$v0<-numeric(length = nrow(velsim))
for(i in 1:nrow(velsim)){
  if(velsim$fac[i]=="1:100"){
    velsim$v0[i] <- coef(nlsMM)[1]*velsim$C_AP[i]/(coef(nlsMM)[4] + velsim$C_AP[i])
  }else{
    if(velsim$fac[i]=="1:200"){
      velsim$v0[i] <- coef(nlsMM)[2]*velsim$C_AP[i]/(coef(nlsMM)[5] + velsim$C_AP[i])
    }else{
      velsim$v0[i] <- coef(nlsMM)[3]*velsim$C_AP[i]/(coef(nlsMM)[6] + velsim$C_AP[i])
    }
  }
}

ggplot(Velocities, aes(C_AP, value))+geom_point(cex=6, pch=21, aes(fill = variable))+
  facet_wrap(~fac, scales="free_y")+theme_min+
  geom_errorbar(aes(ymin = value-se, ymax=value+se))+
  geom_line(data = velsim, aes(C_AP, v0), lwd=1)+
  ylab(expression(paste(italic(v[0]), " (", mu, "mol ", L^{-1}~h^{-1}, ")")))+
  xlab(expression(paste(italic(S[0]), " (", mu, "mol ", L^{-1}, ")")))+
  scale_fill_manual(values = c("white", "grey90", "grey70", "black"),
                    labels = expression(I==0~"& linear range",
                                        I==0~"& full range",
                                        I!=0~"& linear range",
                                        I!=0~"& full range"))+
  theme(legend.title = element_blank(),
        legend.position = c(0.9, 0.25),
        legend.text.align = 0)

#Figure 4: Initial rate of MUB-P substrate decay (v0) at six different initial concentrations of MUB-P 
#(S0). The figure is divided into three boxes representing three different fresh soil - buffer solution 
#slurries (denoted by the box headings). The initial reaction rate was calculated as a slope of increase of a reaction product in 
#time using linear regression with or without the intercept (I = 0 or I 'is not' 0) and time interval spanning 
#either full range or the initial linear interval identified using analysis of residuals 
#(see Material and Methods for details and Fig. S4). Initial reaction rates calculated by linear 
#regression analysis applying different assumptions are denoted by different colors of symbols. 
#The solid black lines represent v0 predicted by non-inhibited Michaelis-Menten kinetic depicted in 
#Fig. 2. Symbols represent mean slope estimates and error bars their standard errors. 
#Note that boxes have different y-axis scales.

#Fig. 5
e5$C_AP2<-round(e5$C_AP, 0)
model5<-nls(Pcorr2~(C_AP-Pblanks)-Km[factor(D)]*lambertW0((C_AP-Pblanks)/Km[factor(D)]*exp((C_AP-Pblanks-Vmax[factor(D)]*Time)/Km[factor(D)])), 
            data = e5[e5$Outlier==0, ],
            start=list(Vmax=c(1/3600, 1/3600, 1/3600), Km=c(10, 10, 10)))
summary(model5)

e5_sim<-data.frame(Time = rep(seq(0, 7*60*60, by=0.1*60*60), times=6),
                   C_AP = rep(rep(unique(e5$C_AP), each=length(seq(0, 7*60*60, by=0.1*60*60))), times = 3),
                   Pblanks = rep(rep(unique(e5$Pblanks), each=length(seq(0, 7*60*60, by=0.1*60*60))), times = 3),
                   D = rep(unique(e5$D), each = length(seq(0, 7*60*60, by=0.1*60*60))*6))
e5_sim$C_AP2<-round(e5_sim$C_AP, 0)
e5_sim$Pred<-predict(model5, newdata = e5_sim)

e5a<-merge(e5[e5$D==1, ], e5d1_twoe$Goodness$Yhat, by=c("Time", "C_AP", "Pcorr2"))
e5b<-merge(e5[e5$D==2, ], e5d2_twoe$Goodness$Yhat, by=c("Time", "C_AP", "Pcorr2"))
e5c<-merge(e5[e5$D==4, ], e5d4_twoe$Goodness$Yhat, by=c("Time", "C_AP", "Pcorr2"))

e5.2<-rbind(e5a, e5b, e5c)

e5.2 %>% group_by(Time, C_AP2, D) %>% 
  summarize(y=mean(Pcorr2, na.rm=T),
            y.sd=sd(Pcorr2, na.rm=T),
            Pred = mean(Pred)) %>%
  ggplot(aes(Time/60/60, y))+geom_point(cex=6, pch=21, aes(fill=factor(D)))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), lwd=1)+
  scale_fill_manual(values=c("white", "grey", "black"), 
                    labels = c("1:100", "1:200", "1:400"))+
  facet_wrap(~C_AP2, scales="free")+theme_min+
  theme(legend.title = element_blank(),
        legend.position = c(0.05, 0.9))+
  xlab("Time (h)")+ylab(expression(paste("Product (", mu, "mol ", L^{-1}, ")")))+
  geom_line(data=e5_sim, aes(x=Time/60/60, y=Pred, color=factor(D)), lwd=1, linetype="dashed")+
  geom_line(aes(x=Time/60/60, y=Pred, color=factor(D)), lwd=1)+
  scale_colour_manual(values = c("grey80", "grey40", "black"), 
                      labels = c("1:100", "1:200", "1:400"))

#Figure 5: Progress curves of MUB-P substrate degradation measured on day 21. 
#The figure is divided into six boxes representing six different initial concentrations (S0) of 
#MUB-P substrate (denoted by the box headings). Empty, grey and black symbols denote three different fresh soil - buffer solution 
#slurries. Dashed lines represent the fit of integrated non-inhibited Michaelis-Menten equation of 
#single enzyme - substrate pair (eq. 4).  Solid lines represent the fit of differential equation 
#with two different enzyme classes - one following non-inhibited Michaelis-Menten kinetic and 
#second following kinetic competitively inhibited by its product. Different line colors denote 
#three different fresh soil - buffer solution slurries, whose data has been fitted by respective 
#equations separately. Symbols represent means and error bars their standard deviations. 
#Note that boxes have different y-axis scales.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Appendix
#Fig. A2
e3$C_AP2<-round(e3$C_AP, 0)
model3<-nls(Pcorr2~(C_AP-Pblanks)-Km[factor(D)]*lambertW0((C_AP-Pblanks)/Km[factor(D)]*exp((C_AP-Pblanks-Vmax[factor(D)]*Time)/Km[factor(D)])), 
            data = e3[e3$Outlier==0, ],
            start=list(Vmax=c(1, 1, 1), Km=c(10, 10, 10)))
summary(model3)
e3_sim<-data.frame(Time = rep(seq(0, 6, by=0.1), times=6),
                   C_AP = rep(rep(unique(e3$C_AP), each=length(seq(0, 6, by=0.1))), times = 3),
                   Pblanks = rep(rep(unique(e3$Pblanks), each=length(seq(0, 6, by=0.1))), times = 3),
                   D = (rep(unique(e3$D), each = length(seq(0, 6, by=0.1))*6)))
e3_sim$C_AP2<-round(e3_sim$C_AP, 0)
e3_sim$Pred<-predict(model3, newdata = e3_sim)

e3 %>% filter(Outlier==0) %>% group_by(Time, C_AP2, D) %>% 
  summarize(y=mean(Pcorr2, na.rm=T),
            y.sd=sd(Pcorr2, na.rm=T)) %>%
  ggplot(aes(Time, y))+geom_point(cex=6, pch=21, aes(fill=factor(D)))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), lwd=1)+
  scale_fill_manual(values=c("white", "grey", "black"), 
                    labels = c("1:100", "1:200", "1:400"))+
  facet_wrap(~C_AP2, scales="free")+theme_min+theme(legend.title = element_blank())+
  xlab("Time (h)")+ylab(expression(paste("Product (", mu, "mol ", L^{-1}, ")")))+
  geom_line(data=e3_sim, aes(x=Time, y=Pred, color=factor(D)), lwd=1)+
  scale_colour_manual(values = c("grey80", "grey40", "black"), 
                      labels = c("1:100", "1:200", "1:400"))+
  theme(legend.title = element_blank(),
        legend.position = c(0.05, 0.9))

#Figure A2: Progress curves of MUB-P substrate degradation measured on day 7. 
#The figure is divided into six boxes representing six different initial concentrations 
#(S0) of MUB-P substrate (denoted by the box headings). Empty, grey and black symbols denote three different fresh 
#soil - buffer solution slurries. Solid lines represent the fit of integrated form of 
#non-inhibited Michaelis-Menten equation (eq. 4 in the main text), which corresponded best to 
#measured data according to non-linear least square regression analysis 
#(see Tab. 2 and section 2.2 in the main text). Different line colors denote three different 
#fresh soil - buffer solution slurries, whose data has been fitted by non-linear regression 
#separately. Symbols represent means and error bars their standard deviations. Note that boxes 
#have different y-axis scales.

#Fig. A3
e3$res<-resid(model3)
e3 %>% filter(Outlier==0) %>% group_by(Time, C_AP2, D) %>% 
  summarize(y=mean(res, na.rm=T),
            y.sd=sd(res, na.rm=T)) %>%
  ggplot(aes(Time, y))+geom_point(cex=6, pch=21, aes(fill=factor(D)))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), lwd=1)+
  scale_fill_manual(values=c("white", "grey", "black"), 
                    labels = c("1:100", "1:200", "1:400"))+
  facet_wrap(~C_AP2, scales="free")+theme_min+
  theme(legend.title = element_blank(),
        legend.position = c(0.05, 0.9))+
  geom_hline(yintercept = 0)+
  xlab("Time (h)")+ylab(expression(paste("Residuals (", mu, "mol ", L^{-1}, ")")))

#Figure A3: Residuals of the fit of integrated form of non-inhibited Michaelis-Menten equation 
#depicted in Fig. A2. The figure is divided into six boxes representing six different initial 
#concentrations (S0) of MUB-P substrate (denoted by the box headings). Empty, grey and black symbols denote three different 
#fresh soil - buffer solution slurries. Solid horizonal line highlight the value 
#zero - i.e. absolute correspondence between the data and the equation fit. 
#Symbols represent means and error bars their standard deviations. Note that boxes have 
#different y-axis scales.

#Fig. A4
e4$fac<-character(length = nrow(e4))
for(i in 1:nrow(e4)){
  if(e4$D[i]==1){
    e4$fac[i]<-c("1:100")
  }else{
    if(e4$D[i]==2){
      e4$fac[i]<-c("1:200")
    }else{
      e4$fac[i]<-c("1:400")
    }
  }
}

ggplot(e4, aes(Time, Pcorr2))+geom_point(cex=6, pch=21)+
  facet_wrap(C_AP2~fac, scales="free")+theme_min+
  geom_abline(aes(intercept = 0, slope=vu), lwd=1, colour="black")+
  geom_abline(aes(intercept = 0, slope=v), lwd=1, colour="black", linetype = "dashed")+
  geom_abline(aes(intercept = Ivu, slope=vuI), lwd=1, colour="grey40")+
  geom_abline(aes(intercept = Iv, slope=vI), lwd=1, colour="grey40", linetype = "dashed")+
  ylab(expression(paste("Product (", mu, "mol ", L^{-1}, ")")))+
  xlab("Time (h)")

#Figure A4: Progress curves of MUB-P substrate degradation measured on day 14. 
#The figure is divided into 12 boxes representing combination of six different initial concentrations 
#(S0, upper labels of the boxes) of MUB-P substrate added to three different fresh soil - buffer 
#solution slurries (lower labels of the boxes). Lines represent fits of different linear regressions.
#Black and grey lines denote linear regression without and without the intercept respectively. 
#Solid lines represent the fit across the whole range of progress curve whereas dashed lines 
#represent the fit across the initial linear region of the progress curve identified using analysis 
#of residuals (for details see Material and methods section). Note that boxes have different y-axis 
#scales. 

#Fig. A5
e6$C_AP2<-round(e6$C_AP, 0)
model6<-nls(Pcorr2~(C_AP-Pblanks)-Km[factor(D)]*lambertW0((C_AP-Pblanks)/Km[factor(D)]*exp((C_AP-Pblanks-Vmax[factor(D)]*Time)/Km[factor(D)])), 
            data = e6[e6$Outlier==0, ],
            start=list(Vmax=c(1/3600, 1/3600, 1/3600), Km=c(10, 10, 10)))
summary(model6)

e6_sim<-data.frame(Time = rep(seq(0, 7*60*60, by=0.1*60*60), times=6),
                   C_AP = rep(rep(unique(e6$C_AP), each=length(seq(0, 7*60*60, by=0.1*60*60))), times = 3),
                   Pblanks = rep(rep(unique(e6$Pblanks), each=length(seq(0, 7*60*60, by=0.1*60*60))), times = 3),
                   D = rep(unique(e6$D), each = length(seq(0, 7*60*60, by=0.1*60*60))*6))
e6_sim$C_AP2<-round(e6_sim$C_AP, 0)
e6_sim$Pred<-predict(model6, newdata = e6_sim)

e6a<-merge(e6[e6$D==1, ], e6d1_twoe$Goodness$Yhat, by=c("Time", "C_AP", "Pcorr2"))
e6b<-merge(e6[e6$D==2, ], e6d2_twoe$Goodness$Yhat, by=c("Time", "C_AP", "Pcorr2"))
e6c<-merge(e6[e6$D==4, ], e6d4_twoe$Goodness$Yhat, by=c("Time", "C_AP", "Pcorr2"))

e6.2<-rbind(e6a, e6b, e6c)

e6.2 %>% group_by(Time, C_AP2, D) %>% 
  summarize(y=mean(Pcorr2, na.rm=T),
            y.sd=sd(Pcorr2, na.rm=T),
            Pred = mean(Pred)) %>%
  ggplot(aes(Time/60/60, y))+geom_point(cex=6, pch=21, aes(fill=factor(D)))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), lwd=1)+
  scale_fill_manual(values=c("white", "grey", "black"), 
                    labels = c("1:100", "1:200", "1:400"))+
  facet_wrap(~C_AP2, scales="free")+theme_min+
  theme(legend.title = element_blank(),
        legend.position = c(0.05, 0.9))+
  xlab("Time (h)")+ylab(expression(paste("Product (", mu, "mol ", L^{-1}, ")")))+
  geom_line(data=e6_sim, aes(x=Time/60/60, y=Pred, color=factor(D)), lwd=1, linetype="dashed")+
  geom_line(aes(x=Time/60/60, y=Pred, color=factor(D)), lwd=1)+
  scale_colour_manual(values = c("grey80", "grey40", "black"), 
                      labels = c("1:100", "1:200", "1:400"))

#Figure A5: Progress curves of MUB-P substrate degradation measured on day 60. 
#The figure is divided into six boxes representing six different initial concentrations (S0) of 
#MUB-P substrate (denoted by the box headings). Empty, grey and black symbols denote three different fresh soil - buffer solution
#slurries. Dashed lines represent the fit of integrated non-inhibited Michaelis-Menten equation of 
#single enzyme - substrate pair (eq. 4 in the main text). Solid lines represent the fit of 
#differential equation with two different enzyme classes - one following non-inhibited 
#Michaelis-Menten kinetic and second following the kinetic competitively inhibited by its product. 
#Different line colors denote three different fresh soil - buffer solution slurries, whose data has 
#been fitted by respective equations separately. Symbols represent means and error bars their 
#standard deviations. Note that boxes have different y-axis scales.

#Fig. A6
##Two enzymes - prediction of the best equation
te_times<-seq(0, 20000, by = 5)

P1<-as.data.frame(ode(twoe_CI, y = c(P = 0, S = unique(e5$C_AP)[1]), 
                      parms = c(Vmax1 = as.numeric(e5d1_twoe$P[1]),
                                Vmax2 = as.numeric(e5d1_twoe$P[2]),
                                Km = as.numeric(e5d1_twoe$P[3]),
                                Kic = as.numeric(e5d1_twoe$P[4])), times = te_times))
P1$S<-unique(e5$C_AP)[1]

for(i in unique(e5$C_AP)[-1]){
  Px <- as.data.frame(ode(twoe_CI, y = c(P = 0, S = i), 
                                    parms = c(Vmax1 = as.numeric(e5d1_twoe$P[1]),
                                              Vmax2 = as.numeric(e5d1_twoe$P[2]),
                                              Km = as.numeric(e5d1_twoe$P[3]),
                                              Kic = as.numeric(e5d1_twoe$P[4])), times = te_times))
  Px$S <- i
  P1 <- rbind(P1, Px)
}

P1$Legend <- c("No substrate inhibition")

##Two enzymes - prediction of the best equation with substrate inhibition term 
##identified during the first enzyme assay
twoe_CI2<-function(time, state, pars){
  
  with(as.list(c(state, pars)),{
    
    dP<-(Vmax1*S/(Km1 + S + S^2/Ks))+(Vmax2*S/(Km2*(1+P/Kic) + S))
    dS<--(Vmax1*S/(Km1 + S + S^2/Ks))-(Vmax2*S/(Km2*(1+P/Kic) + S))
    
    return(list(c(dP, dS)))
    
  })
}

P2<-as.data.frame(ode(twoe_CI2, y = c(P = 0, S = unique(e5$C_AP)[1]), 
                      parms = c(Vmax1 = as.numeric(e5d1_twoe$P[1]),
                                Vmax2 = as.numeric(e5d1_twoe$P[2]),
                                Km2 = as.numeric(e5d1_twoe$P[3]),
                                Kic = as.numeric(e5d1_twoe$P[4]),
                                Km1 = 17.92652, Ks = 215.5652), times = te_times))
P2$S<-unique(e5$C_AP)[1]

for(i in unique(e5$C_AP)[-1]){
  Px <- as.data.frame(ode(twoe_CI2, y = c(P = 0, S = i), 
                          parms = c(Vmax1 = as.numeric(e5d1_twoe$P[1]),
                                    Vmax2 = as.numeric(e5d1_twoe$P[2]),
                                    Km2 = as.numeric(e5d1_twoe$P[3]),
                                    Kic = as.numeric(e5d1_twoe$P[4]),
                                    Km1 = 17.92652, Ks = 215.5652), times = te_times))
  Px$S <- i
  P2 <- rbind(P2, Px)
}

P2$Legend <- c("Substrate inhibition")

ggplot(Ps, aes(time/60/60, P))+
  geom_line(aes(color=factor(S), linetype = Legend), lwd = 1.5, alpha=0.5)+
  theme_min+xlab("Time (h)")+ylab(expression(paste("Product (", mu, "mol ", L^{-1}, ")")))+
  theme(legend.title = element_blank(),
        legend.text.align = 0)+
  scale_color_manual(values = c("black", "grey30", "grey80", "dodgerblue4", "dodgerblue2",
                                "cornflowerblue"),
                        labels = expression(S==2~mu*mol~L^{-1},
                                            S==10~mu*mol~L^{-1},
                                            S==21~mu*mol~L^{-1}, 
                                            S==52~mu*mol~L^{-1},
                                            S==104~mu*mol~L^{-1},
                                            S==207~mu*mol~L^{-1}))

#Figure A6: Simulation of the progress curves for six different initial concentrations of 
#MUB-P substrate (S0) measured at 1:100 soil - buffer slurry. Different S0 are denoted by 
#different line colours. Solid lines represent the simulation of the differential equation 5, 
#which assumes formation of the product by combined activity of two acid phosphatases pools, 
#one of which being inhibited by the product. The parameters of eq. 5 used for this simulation 
#are reported in Table A2 (i.e. for day 21). The dashed lines represent the fit of the same 
#equation but assuming the substrate inhibition of the second acid phosphatase pool with 
#inhibition constant Ki = 215.6 estimated for the acid phosphatase activity on day zero.

###########################################For reviewers################################
fdata<-data.frame(Time=c(0.05, 1, 0, 2),
                  P= c(0, 6, 0, 8),
                  Group=c("A", "A", "B", "B"))

ggplot(fdata, aes(Time, P)) + geom_point(size=8, pch=21, aes(fill = Group),
                                         show.legend = F)+
  theme_min+theme(axis.text.y = element_blank())+
  ylab("Product concentration")+
  xlab("Time (hours)")+
  stat_smooth(method = lm, aes(color=Group), show.legend = F)+
  scale_x_continuous(limits=c(0,2), breaks=c(0,1,2))
  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
############################################END OF SCRIPT#####################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
