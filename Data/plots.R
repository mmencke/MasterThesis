#########################################
# BASIC SETUP                           #
#########################################

rm(list=ls())
library(ggplot2)
library(extrafont)
library(gridExtra)
library(reshape2)
#library(gg3D)

#Set WD to current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


folder<-"CppOutput"

load.matrix<-function(folder_name,file_name) {
  file_path<-paste(folder_name,file_name,sep="/")
  return(data.table::fread(file_path,header=F,sep=","))
}


#########################################
# COLOURS                               #
#########################################


#primary colours
cbs.blue=rgb(73,103,170,max=255)
cbs.black=rgb(51,51,51,max=255)
cbs.darkgrey=rgb(112,112,112,max=255)
cbs.background=rgb(247,247,243,max=255)

#secondary colours
cbs.yellow=rgb(253,207,65,max=255)
cbs.orange=rgb(255,150,100,max=255)
cbs.red=rgb(255,102,94,max=255)
cbs.purple=rgb(133,134,198,max=255)

cbs.turquoise=rgb(72,184,231,max=255)
cbs.green=rgb(80,214,145,max=255)
cbs.brown=rgb(210,162,124,max=255)
cbs.lightgrey=rgb(185,194,200,max=255)

theme_cbs <- function(){ 
  
  theme_minimal() %+replace%    #replace elements we want to change
    
    theme(
      panel.background = element_rect(fill = cbs.background,
                                      colour = cbs.lightgrey,
                                      size = 0.5, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      colour = cbs.lightgrey), 
      panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                      colour = cbs.lightgrey),
      text=element_text(size=20, family="Adobe Garamond Pro")
    )
}


palette(c(cbs.blue,cbs.yellow,cbs.red,cbs.purple,cbs.turquoise,cbs.green))

#########################################
# CORRELATED NORMALS                    #
#########################################

correlated_normals<-load.matrix(folder, "correlatedNormals.csv")


p_cor_norm_ind<-ggplot(data=correlated_normals) +  
  geom_point(aes(x=V1,y=V2),colour=cbs.blue)+
  labs(x="Variable 1", y="Variable 2",title="Uncorrelated")+
  theme_cbs()

p_cor_norm_pos<-ggplot(data=correlated_normals) +  
  geom_point(aes(x=V1,y=V3),colour=cbs.blue)+
  labs(x="Variable 1", y="Variable 3",title="Positive Correlation")+
  theme_cbs()

p_cor_norm_neg<-ggplot(data=correlated_normals) +  
  geom_point(aes(x=V2,y=V3),colour=cbs.blue)+
  labs(x="Variable 2", y="Variable 3",title="Negative Correlation")+
  theme_cbs()


ggsave("ROutput/correlated-normals.pdf", arrangeGrob(p_cor_norm_ind, p_cor_norm_pos, p_cor_norm_neg, nrow=2),width=16, height = 9)
embed_fonts("ROutput/correlated-normals.pdf")



#########################################
# VASICEK DISTRIBUTION                  #
#########################################

vasicek_dist<-load.matrix(folder, "vasicekDistribution.csv")

kappa<-0.05
theta<-0.03
sigma<-0.01
r0<-0.03
t<-1

#From [Andersen and Piterbarg (2010), p. 413]
vasicek_mean <- theta+(r0-theta)*exp(-kappa*t)
vasicek_sd<-sqrt(sigma^2/(2*kappa)*(1-exp(-2*kappa*t)))

xfit <- seq(min(vasicek_dist$V1), max(vasicek_dist$V1), length = 10000) 
yfit <- dnorm(xfit, mean = vasicek_mean, sd = vasicek_sd) 

true_density<-data.frame(xfit,yfit)


p_vasicek<-ggplot() + 
  geom_density(data=vasicek_dist,
               aes(V1,
                   colour="Discretised"),
               size=1.5,linetype=1,show.legend = F)+ 
  geom_line(data=true_density,
            aes(x=xfit,y=yfit,
                colour="True",),
            size=1.5) +
  labs(colour="Legend", x="Interest Rate", y="Density")+
  scale_colour_manual(values=c("Discretised"=cbs.black, "True"=cbs.blue))+
  theme_cbs()


ggsave("ROutput/vasicek-dist.pdf", plot=p_vasicek,  width=16, height=9)

embed_fonts("ROutput/vasicek-dist.pdf")

#########################################
# COMPARE CIR DISTRIBUTION              #
#########################################

cir_dist<-load.matrix(folder, "cirDistribution.csv")

t<-1
kappa_cir<-0.4
theta_cir<-0.026
sigma_cir<-0.14
r0_cir<-0.0165

#From [Brigo, Morini and Pallavicini (2013), p. 73]
c<- 4*kappa_cir/(sigma_cir^2*(1-exp(-kappa_cir*t)))
v<- 4*kappa_cir*theta_cir/sigma_cir^2
lambda<- c*r0_cir*exp(-kappa_cir*t)

xfit <- seq(min(cir_dist$V5)-0.01, max(cir_dist$V5)+0.01, length = 10000)
yfit <- c*dchisq(xfit*c,v,ncp=lambda) 
true_density<-data.frame(xfit,yfit)

p_cir_simple<-ggplot() + 
  geom_density(data=cir_dist,
               aes(V1,
                   colour="Discretised"),
               size=1.5,linetype=1,show.legend = F)+ 
  geom_line(data=true_density,
            aes(x=xfit,y=yfit,
                colour="True",),
            size=1.5) +
  labs(colour="Legend", x="Interest Rate", y="Density",title="Simple Euler")+
  scale_colour_manual(values=c("Discretised"=cbs.black, "True"=cbs.blue))+
  theme_cbs()

p_cir_trunc<-ggplot() + 
  geom_density(data=cir_dist,
               aes(V2,
                   colour="Discretised"),
               size=1.5,linetype=1,show.legend = F)+ 
  geom_line(data=true_density,
            aes(x=xfit,y=yfit,
                colour="True",),
            size=1.5) +
  labs(colour="Legend", x="Interest Rate", y="Density",title="Full Truncation")+
  scale_colour_manual(values=c("Discretised"=cbs.black, "True"=cbs.blue))+
  theme_cbs()


p_cir_sqrt<-ggplot() + 
  geom_density(data=cir_dist,
               aes(V3,
                   colour="Discretised"),
               size=1.5,linetype=1,show.legend = F)+ 
  geom_line(data=true_density,
            aes(x=xfit,y=yfit,
                colour="True",),
            size=1.5) +
  labs(colour="Legend", x="Interest Rate", y="Density",title="Simulate Square Root of CIR")+
  scale_colour_manual(values=c("Discretised"=cbs.black, "True"=cbs.blue))+
  theme_cbs()


p_cir_qe<-ggplot() + 
  geom_density(data=cir_dist,
               aes(V4,
                   colour="Discretised"),
               size=1.5,linetype=1,show.legend = F)+ 
  geom_line(data=true_density,
            aes(x=xfit,y=yfit,
                colour="True",),
            size=1.5) +
  labs(colour="Legend", x="Interest Rate", y="Density",title="Quadratic Exponential")+
  scale_colour_manual(values=c("Discretised"=cbs.black, "True"=cbs.blue))+
  theme_cbs()

p_cir_exact<-ggplot() + 
  geom_density(data=cir_dist,
               aes(V5,
                   colour="Discretised"),
               size=1.5,linetype=1,show.legend = F)+ 
  geom_line(data=true_density,
            aes(x=xfit,y=yfit,
                colour="True",),
            size=1.5) +
  labs(colour="Legend", x="Interest Rate", y="Density",title="Exact")+
  scale_colour_manual(values=c("Discretised"=cbs.black, "True"=cbs.blue))+
  theme_cbs()


ggsave("ROutput/cir-dist.pdf", arrangeGrob(p_cir_simple, p_cir_trunc,p_cir_sqrt,p_cir_qe,p_cir_exact),width=16, height = 9)
embed_fonts("ROutput/cir-dist.pdf")

cir_qe_convergence<-load.matrix(folder, "cirQeConvergence.csv")

xfit <- seq(min(cir_qe_convergence$V1), max(cir_qe_convergence$V1), length = 1000000)
yfit <- c*dchisq(xfit*c,v,ncp=lambda) 
true_density<-data.frame(xfit,yfit)


p_cir_qe_conv<-ggplot() + 
  geom_density(data=cir_qe_convergence,
               aes(V1,
                   colour="Discretised"),
               size=1.5,linetype=1,show.legend = F)+ 
  geom_line(data=true_density,
            aes(x=xfit,y=yfit,
                colour="True",),
            size=1.5) +
  labs(colour="Legend", x="Interest Rate", y="Density",title="Quadratic Exponential")+
  scale_colour_manual(values=c("Discretised"=cbs.black, "True"=cbs.blue))+
  theme_cbs()


ggsave("ROutput/cir-qe-convergence.pdf", p_cir_qe_conv, width=16, height = 9)
embed_fonts("ROutput/cir-qe-convergence.pdf")

#########################################
# CALIBRATION OF CIR TO OLD MARKET DATA #
#########################################

calibrated_cir_term<-load.matrix(folder, "calibratedCirTerm.csv")
p_cir_cal_cds<-ggplot(data=calibrated_cir_term) +  
  geom_point(aes(x=V1,y=V3,colour="CIR"))+
  geom_point(aes(x=V1,y=V4,colour="CIR++"))+
  geom_line(aes(x=V1,y=V2))+
  scale_colour_manual(values=c("True"=cbs.black,"CIR"=cbs.blue,"CIR++"=cbs.yellow))+
  labs(colour="Legend",x="Time", y="Survival Probability")+
  theme_cbs()

calibrated_cir<-load.matrix(folder, "calibratedCir.csv")
p_cir_cal_cdso<-ggplot(data=calibrated_cir) +  
  geom_point(aes(x=V3,y=V4,colour="CIR"))+
  geom_point(aes(x=V3,y=V5,colour="CIR++"))+
  geom_abline(aes(x=V3),slope=1,intercept=0,colour=cbs.black)+
  scale_colour_manual(values=c("CIR"=cbs.blue,"CIR++"=cbs.yellow))+
  labs(colour="Legend",x="Market Price", y="Model Price")+
  theme_cbs()

ggsave("ROutput/calibrated-cir.pdf", arrangeGrob(p_cir_cal_cds, p_cir_cal_cdso,nrow=1),width=16, height = 9)
embed_fonts("ROutput/calibrated-cir.pdf")

#########################################
# SWAP RATES                            #
#########################################

swap_tenors<-c(1,2,3,4,5,6,7,8,9,10,11,12,15,20,25,30,40,50)
swap_rates<-c(-0.5425,-0.5224,-0.479,-0.413,-0.3407,-0.2605,-0.179,-0.099,-0.022,0.0505,0.1175,0.1795,0.333,0.4645,0.4994,0.494,0.459,0.4227)

swap_data<-data.frame(swap_tenors,swap_rates)
p_swap<-ggplot(data=swap_data,aes(x=swap_tenors,y=swap_rates)) +
  geom_point(colour=cbs.blue,size=3)+theme_cbs()+
  labs(x="Swap Tenor", y="Swap Rate")

ggsave("ROutput/swap-rates.pdf",p_swap,width=16, height = 9)
embed_fonts("ROutput/swap-rates.pdf")



#########################################
# CALIBRATED ZERO CURVE                 #
#########################################

calibratedZeroCurve<-load.matrix(folder, "calibratedZeroCurve.csv")
p_cal_zero<-ggplot(data=calibratedZeroCurve) +  
  geom_point(aes(x=V1,y=V3,colour="G2++"))+
  geom_point(aes(x=V1,y=V4,colour="Hull-White"))+
  geom_line(aes(x=V1,y=V2))+
  scale_colour_manual(values=c("G2++"=cbs.blue, "Hull-White"=cbs.yellow))+
  labs(colour="Legend",x="Time", y="Zero Rate")+
  theme_cbs()

ggsave("ROutput/calibrated-zero-curve.pdf", p_cal_zero,width=16, height = 9)
embed_fonts("ROutput/calibrated-zero-curve.pdf")



#########################################
# SWAPTION VOLS                         #
#########################################


swaption_vols<-load.matrix("CppInput", "swaptionVols.csv")

swaption_expiries<-c("1Mo", "2Mo", "3Mo", "6Mo", "9Mo", "1Yr", "18Mo", "2Yr", "3Yr", "4Yr", "5Yr", "6Yr", "7Yr", "8Yr", "9Yr", "10Yr", "15Yr", "20Yr", "25Yr", "30Yr")
swaption_tenors<-c("1Yr", "2Yr", "3Yr", "4Yr", "5Yr", "6Yr", "7Yr", "8Yr", "9Yr", "10Yr", "15Yr", "20Yr", "25Yr", "30Yr")


swaption_vols_melt<-data.frame()
k<-1
for(i in 1:dim(swaption_vols)[1]) {
  for(j in 1:dim(swaption_vols)[2]) {
    swaption_vols_melt[k,1]<-swaption_expiries[i]
    swaption_vols_melt[k,2]<-swaption_tenors[j]
    swaption_vols_melt[k,3]<-round(swaption_vols[[j]][[i]],0)
      
    k<-k+1
  }
}

p_swaption<-ggplot(data=swaption_vols_melt,aes(V1,V2)) +
  geom_tile(aes(fill = V3)) +
  geom_text(aes(label=V3),family="Adobe Garamond Pro")+
  scale_x_discrete(limits=swaption_expiries)+
  scale_y_discrete(limits=swaption_tenors)+
  scale_fill_gradient2(mid=cbs.yellow,
                       high = cbs.orange)+
  theme_cbs()+
  labs(x="Expiry",y="Tenor",fill="Volatility")

ggsave("ROutput/swaption-vol-surface.pdf",p_swaption,width=16, height = 9)
embed_fonts("ROutput/swaption-vol-surface.pdf")

#########################################
# SWAPTION FIT                          #
#########################################
swaption_fit<-load.matrix(folder, "swaptionFit.csv")

swaption_fit_round<-swaption_fit
swaption_fit_round$V1<-factor(as.character(round(swaption_fit_round$V1,2)))
swaption_fit_round$V2<-factor(as.character(round(swaption_fit_round$V2,2)))

swaption_fit_round$V6<-abs(swaption_fit_round$V4-swaption_fit_round$V3)#G2
swaption_fit_round$V7<-abs(swaption_fit_round$V5-swaption_fit_round$V3)#Hull-White

sum(swaption_fit_round$V6)#G2
sum(swaption_fit_round$V7)#Hull-White

swaption_fit_round$V3<-round(swaption_fit_round$V3)
swaption_fit_round$V4<-round(swaption_fit_round$V4)
swaption_fit_round$V5<-round(swaption_fit_round$V5)
swaption_fit_round$V6<-round(swaption_fit_round$V6)
swaption_fit_round$V7<-round(swaption_fit_round$V7)

x_lim<-unique(swaption_fit_round$V1)
y_lim<-unique(swaption_fit_round$V2)

p_swaption_market<-ggplot(data=swaption_fit_round,aes(V1,V2)) +
  geom_tile(aes(fill = V3)) +
  geom_text(aes(label=V3),family="Adobe Garamond Pro")+
  scale_x_discrete(limits=x_lim)+
  scale_y_discrete(limits=y_lim)+
  scale_fill_gradient2(mid=cbs.yellow,
                       high = cbs.orange)+
  theme_cbs()+
  labs(x="Expiry",y="Tenor",fill="Volatility")


p_swaption_g2<-ggplot(data=swaption_fit_round,aes(V1,V2)) +
  geom_tile(aes(fill = V4)) +
  geom_text(aes(label=V4),family="Adobe Garamond Pro")+
  scale_x_discrete(limits=x_lim)+
  scale_y_discrete(limits=y_lim)+
  scale_fill_gradient2(mid=cbs.yellow,
                       high = cbs.orange)+
  theme_cbs()+
  labs(x="Expiry",y="Tenor",fill="Volatility")

p_swaption_g2_dif<-ggplot(data=swaption_fit_round,aes(V1,V2)) +
  geom_tile(aes(fill = V6)) +
  geom_text(aes(label=V6),family="Adobe Garamond Pro")+
  scale_x_discrete(limits=x_lim)+
  scale_y_discrete(limits=y_lim)+
  scale_fill_gradient2(mid=cbs.yellow,
                       high = cbs.orange)+
  theme_cbs()+
  labs(x="Expiry",y="Tenor",fill="Abs. Dif.")

ggsave("ROutput/g2-calibration.pdf",arrangeGrob(p_swaption_g2, p_swaption_g2_dif,nrow=1),width=22, height = 9)
embed_fonts("ROutput/g2-calibration.pdf")

p_swaption_hw<-ggplot(data=swaption_fit_round,aes(V1,V2)) +
  geom_tile(aes(fill = V5)) +
  geom_text(aes(label=V5),family="Adobe Garamond Pro")+
  scale_x_discrete(limits=x_lim)+
  scale_y_discrete(limits=y_lim)+
  scale_fill_gradient2(mid=cbs.yellow,
                       high = cbs.orange)+
  theme_cbs()+
  labs(x="Expiry",y="Tenor",fill="Volatility")

p_swaption_hw_dif<-ggplot(data=swaption_fit_round,aes(V1,V2)) +
  geom_tile(aes(fill = V7)) +
  geom_text(aes(label=V7),family="Adobe Garamond Pro")+
  scale_x_discrete(limits=x_lim)+
  scale_y_discrete(limits=y_lim)+
  scale_fill_gradient2(mid=cbs.yellow,
                       high = cbs.orange)+
  theme_cbs()+
  labs(x="Expiry",y="Tenor",fill="Abs. Dif.")


ggsave("ROutput/hw-calibration.pdf",arrangeGrob(p_swaption_hw, p_swaption_hw_dif,nrow=1),width=22, height = 9)
embed_fonts("ROutput/hw-calibration.pdf")



#########################################
# CREDIT SECTOR CURVES                  #
#########################################



cdsSectorCurves<-data.table::fread("RInput/cdsSectorCurves.csv",header=T,sep=",")

financials<-cdsSectorCurves[c(cdsSectorCurves$Region=="Europe"),]
financials<-financials[financials$Sector=="Financials",]
financials<-financials[,c(5,8:18)]
colnames(financials)<-c("Rating",0.5,1:5,7,10,15,20,30)
financials<-melt(financials)

p_financials<-ggplot(data=financials,aes(x=variable,y=value,colour=Rating))+
  geom_point(size=3)+
  theme_cbs()+
  labs(x="Maturity",y="Spread",title="Financials")+
  scale_colour_manual(values=c(cbs.blue,cbs.yellow,cbs.red,cbs.purple,cbs.green,cbs.turquoise,cbs.orange))

corporates<-cdsSectorCurves[c(cdsSectorCurves$Region=="Europe"),]
corporates<-corporates[corporates$Sector=="Corporates",]
corporates<-corporates[,c(5,8:18)]
colnames(corporates)<-c("Rating",0.5,1:5,7,10,15,20,30)
corporates<-melt(corporates)

p_corporates<-ggplot(data=corporates,aes(x=variable,y=value,colour=Rating))+
  geom_point(size=3)+
  theme_cbs()+
  labs(x="Maturity",y="Spread",title="Corporates")+
  scale_colour_manual(values=c(cbs.blue,cbs.yellow,cbs.red,cbs.purple,cbs.green,cbs.turquoise,cbs.orange))

ggsave("ROutput/cds-sector-curves.pdf",arrangeGrob(p_financials, p_corporates,nrow=1),width=16, height = 9)
embed_fonts("ROutput/cds-sector-curves.pdf")



#########################################
# CALIBRATED DEFAULT CURVE              #
#########################################


calibrated_default_curve<-load.matrix(folder, "calibratedDefaultCurve.csv")
p_cal_default_bank<-ggplot(data=calibrated_default_curve) +  
  geom_point(aes(x=V1,y=V3,colour="CIR++"))+
  geom_line(aes(x=V1,y=V2))+
  scale_colour_manual(values=c("CIR++"=cbs.blue))+
  labs(colour="Legend",x="Time", y="Surv. Prob.",title="Bank")+
  theme_cbs()
p_cal_default_cpt<-ggplot(data=calibrated_default_curve) +  
  geom_point(aes(x=V1,y=V5,colour="CIR++"))+
  geom_line(aes(x=V1,y=V4))+
  scale_colour_manual(values=c("CIR++"=cbs.blue))+
  labs(colour="Legend",x="Time", y="Surv. Prob.",title="Counterparty")+
  theme_cbs()

ggsave("ROutput/calibrated-default-curve.pdf", arrangeGrob(p_cal_default_bank,p_cal_default_cpt,nrow=1),width=16, height = 9)
embed_fonts("ROutput/calibrated-default-curve.pdf")


#########################################
# CREDIT SPREAD VOLATILITY              #
#########################################


cds_vol_dates<-c(43861,43889,43921,43951,43980,44012,44043,44071,44104,44134,44165,44196,44225,44253,44286,44316,44321)
cds_vol_dates<-as.Date(cds_vol_dates, origin = "1899-12-30")
cds_vols<-c(0.29553358150798,0.427535610455252,0.960886456441857,1.0769910207602,1.10271650757223,1.14960761931817,0.779261020200403,0.505797414663751,0.496654847484805,0.446810728951232,0.462661932100958,0.503450911373786,0.479622113351782,0.452303369069778,0.378740011411154,0.332739537632943,0.305996550075626)

cds_vol<-data.frame(cds_vol_dates,cds_vols)

p_cds_vol<-ggplot(data=cds_vol,aes(x=cds_vol_dates,y=cds_vols)) +
  geom_point(colour=cbs.blue,size=3)+theme_cbs()+
  labs(x="Date", y="Volatility")

ggsave("ROutput/cds-vol.pdf",p_cds_vol,width=16, height = 9)
embed_fonts("ROutput/cds-vol.pdf")

#########################################
# HOSTRUPS HAVE                         #
#########################################

tenor<-c(1,2,3,4,5,6,7,8,9,10,11,12,15,20,25,30,40,50)
rates2006<-c(4.07,4.12,4.1263,4.126,4.1225,4.1305,4.1426,4.1586,4.1771,4.1961,4.2131,4.2301,4.2724,4.3084,4.3064,4.287,4.244,4.202)
rates2016<-c(-0.302,-0.2695,-0.2195,-0.148,-0.0535,0.06,0.1875,0.315,0.438,0.543,0.653,0.7447,0.948,1.1149,1.156,1.194,1.186,1.134)

hostrupHave<-data.frame(tenor,rates2006,rates2016)

p_hostrupHave<-ggplot(data=hostrupHave) +  
  geom_point(aes(x=tenor,y=rates2006,colour="29 Dec 2006"),size=3)+
  geom_point(aes(x=tenor,y=rates2016,colour="30 Dec 2016"),size=3)+
  scale_colour_manual(values=c("29 Dec 2006"=cbs.blue,"30 Dec 2016"=cbs.yellow))+
  labs(colour="Reference Date",x="Maturity", y="Swap Rate",title="Term Structure")+
  theme_cbs()


hhDates<-as.Date(c(39080,39447,39813,40178,40543,40907,41274,41639,42004,42369,42734),origin = "1899-12-30")
hh20y<-c(4.3084,4.9,3.812,3.972,3.615,2.622,2.073,2.64,1.239,1.5073,1.1149)
hh15y<-c(4.2724,4.848,3.8465,3.865,3.54,2.5835,1.903,2.492,1.047,1.3231,0.948)
hh12y<-c(4.2301,4.774,3.782,3.675,3.354,2.44,1.668,2.264,0.86,1.103,0.7447)

hhDevelopment<-data.frame(hhDates,hh12y,hh15y,hh20y)

p_hhDevelopment<-ggplot(data=hhDevelopment) +  
  geom_point(aes(x=hhDates,y=hh12y,colour="12y"),size=3)+
  geom_point(aes(x=hhDates,y=hh15y,colour="15y"),size=3)+
  geom_point(aes(x=hhDates,y=hh20y,colour="20y"),size=3)+
  scale_colour_manual(values=c("12y"=cbs.blue,"15y"=cbs.yellow,"20y"=cbs.red))+
  labs(colour="Maturity",x="Reference Date", y="Swap Rate",title="Development")+
  theme_cbs()

ggsave("hostrup-have.pdf", 
       arrangeGrob(p_hostrupHave, p_hhDevelopment, nrow=1),
       width=16, height = 9)
embed_fonts("hostrup-have.pdf")

#########################################
# VARIOUS CALCULATIONS                  #
#########################################

#VARIANCE IN CIR
cirVar<-function(y0,kappa,theta,sigma,t) {
  exp.kappa<-exp(-kappa*t)
  
  var<-(y0*sigma^2*exp.kappa)/kappa*(1-exp.kappa)+(theta*sigma^2)/(2*kappa)*(1-exp.kappa)^2
return(var)
  }

sqrt(cirVar(0.0165,0.4,0.026,0.14,1))



#95% CONFIDENCE INTERVAL (2-SIDED)
sprintf("%.15f",qnorm(0.975))


#FIND CORRELATION BETWEEN INTENSITY AND G2
sigma1<- 0.00970193
sigma2<- 0.00955238
rho12<- -0.759051

std.dev<-sqrt(sigma1^2+sigma2^2+2*sigma1*sigma2*rho12)
f<-function(rho13,rho23) {
  return((sigma1*rho13+sigma2*rho23)/std.dev)
}

fToSolve<-function(rho23,c) {
  return(f(-0.05,rho23)-c)
}

uniroot(fToSolve, c(-1,1),c=1)$root
uniroot(fToSolve, c(-1,1),c=-1)$root

sprintf("%.15f",f(-0.05,-0.648994))
sprintf("%.15f",f(-0.05,0.7505596))


#########################################
# TRANSFORM DATA TO TABLES              #
#########################################

exmp_swaps<-load.matrix(folder, "exmpSwaps.csv")

exmp_swaps_output<-data.frame(matrix(NA,nrow=dim(exmp_swaps)[1],ncol=dim(exmp_swaps)[2]))

for(i in 1:dim(exmp_swaps)[1]) {
  exmp_swaps_output$X1[i]<-exmp_swaps$V1[i]
  exmp_swaps_output$X2[i]<-format(round(exmp_swaps$V2[i]),big.mark = ",")
  exmp_swaps_output$X3[i]<-format(round(exmp_swaps$V3[i]),big.mark = ",")
  exmp_swaps_output$X4[i]<-format(round(exmp_swaps$V4[i]),big.mark = ",")
  exmp_swaps_output$X5[i]<-format(round(exmp_swaps$V5[i]),big.mark = ",")
  exmp_swaps_output$X6[i]<-format(round(exmp_swaps$V6[i]),big.mark = ",")
  exmp_swaps_output$X7[i]<-format(round(exmp_swaps$V7[i]),big.mark = ",")
}
write.table(exmp_swaps_output, file="ROutput/exmp-swaps-output.csv",
            quote=F, sep=";", col.names = F,row.names = F)


cva_dva_ind<-load.matrix(folder, "cvaDvaInd.csv")

cva_dva_ind_output<-data.frame(matrix(NA,nrow=6,ncol=7))

for(i in 1:dim(cva_dva_ind)[1]) {
  cva_dva_ind_output$X1[i]<-cva_dva_ind$V1[i]
  cva_dva_ind_output$X2[i]<-paste0(format(round(cva_dva_ind$V2[i]),big.mark = ",")," (", format(round(cva_dva_ind$V3[i]),big.mark=","),")")
  cva_dva_ind_output$X3[i]<-paste0(format(round(cva_dva_ind$V4[i]),big.mark=",")," (", format(round(cva_dva_ind$V5[i]),big.mark=","),")")
  cva_dva_ind_output$X4[i]<-paste0(format(round(cva_dva_ind$V6[i]),big.mark=",")," (", format(round(cva_dva_ind$V7[i]),big.mark=","),")")
  cva_dva_ind_output$X5[i]<-paste0(format(round(cva_dva_ind$V8[i]),big.mark=",")," (", format(round(cva_dva_ind$V9[i]),big.mark=","),")")
  cva_dva_ind_output$X6[i]<-paste0(format(round(cva_dva_ind$V10[i]),big.mark=",")," (", format(round(cva_dva_ind$V11[i]),big.mark=","),")")
  cva_dva_ind_output$X7[i]<-paste0(format(round(cva_dva_ind$V12[i]),big.mark=",")," (", format(round(cva_dva_ind$V13[i]),big.mark=","),")")
}

write.table(cva_dva_ind_output, file="ROutput/cva-dva-ind-output.csv",
            quote=F, sep=";", col.names = F,row.names = F)


cva_dva_wwr<-load.matrix(folder, "cvaDvaWwr.csv")

cva_dva_wwr_output<-data.frame(matrix(NA,nrow=6,ncol=7))

for(i in 1:dim(cva_dva_wwr)[1]) {
  cva_dva_wwr_output$X1[i]<-cva_dva_wwr$V1[i]
  cva_dva_wwr_output$X2[i]<-paste0(format(round(cva_dva_wwr$V2[i]),big.mark = ",")," (", format(round(cva_dva_wwr$V3[i]),big.mark=","),")")
  cva_dva_wwr_output$X3[i]<-paste0(format(round(cva_dva_wwr$V4[i]),big.mark=",")," (", format(round(cva_dva_wwr$V5[i]),big.mark=","),")")
  cva_dva_wwr_output$X4[i]<-paste0(format(round(cva_dva_wwr$V6[i]),big.mark=",")," (", format(round(cva_dva_wwr$V7[i]),big.mark=","),")")
  cva_dva_wwr_output$X5[i]<-paste0(format(round(cva_dva_wwr$V8[i]),big.mark=",")," (", format(round(cva_dva_wwr$V9[i]),big.mark=","),")")
  cva_dva_wwr_output$X6[i]<-paste0(format(round(cva_dva_wwr$V10[i]),big.mark=",")," (", format(round(cva_dva_wwr$V11[i]),big.mark=","),")")
  cva_dva_wwr_output$X7[i]<-paste0(format(round(cva_dva_wwr$V12[i]),big.mark=",")," (", format(round(cva_dva_wwr$V13[i]),big.mark=","),")")
}

write.table(cva_dva_wwr_output, file="ROutput/cva-dva-wwr-output.csv",
            quote=F, sep=";", col.names = F,row.names = F)

cva_dva_rwr<-load.matrix(folder, "cvaDvaRwr.csv")

cva_dva_rwr_output<-data.frame(matrix(NA,nrow=6,ncol=7))

for(i in 1:dim(cva_dva_rwr)[1]) {
  cva_dva_rwr_output$X1[i]<-cva_dva_rwr$V1[i]
  cva_dva_rwr_output$X2[i]<-paste0(format(round(cva_dva_rwr$V2[i]),big.mark = ",")," (", format(round(cva_dva_rwr$V3[i]),big.mark=","),")")
  cva_dva_rwr_output$X3[i]<-paste0(format(round(cva_dva_rwr$V4[i]),big.mark=",")," (", format(round(cva_dva_rwr$V5[i]),big.mark=","),")")
  cva_dva_rwr_output$X4[i]<-paste0(format(round(cva_dva_rwr$V6[i]),big.mark=",")," (", format(round(cva_dva_rwr$V7[i]),big.mark=","),")")
  cva_dva_rwr_output$X5[i]<-paste0(format(round(cva_dva_rwr$V8[i]),big.mark=",")," (", format(round(cva_dva_rwr$V9[i]),big.mark=","),")")
  cva_dva_rwr_output$X6[i]<-paste0(format(round(cva_dva_rwr$V10[i]),big.mark=",")," (", format(round(cva_dva_rwr$V11[i]),big.mark=","),")")
  cva_dva_rwr_output$X7[i]<-paste0(format(round(cva_dva_rwr$V12[i]),big.mark=",")," (", format(round(cva_dva_rwr$V13[i]),big.mark=","),")")
}

write.table(cva_dva_rwr_output, file="ROutput/cva_dva_rwr_output.csv",
            quote=F, sep=";", col.names = F,row.names = F)

