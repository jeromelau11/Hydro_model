
#################################################################################
# Model information
# Model name: HBV-Extended (HBV-E) hydrological model
# Model developer: Zhiyong Liu 
# Contact: liuzhiy25@mail.sysu.edu.cn
# Version: HBV-E v1.1
# Language: R 
# Model type: physically-based / conceptual 
# Time scales: applications in daily time steps (aims to hourly and monthly scales)
# Update Date: 2025-04-01
#################################################################################


library(xts)
library(hydroTSM)
library(hydroGOF)
library(Evapotranspiration)
require(lubridate)
library(dream)
library(gsubfn)

orginal_path<-"github\\data\\"
setwd(orginal_path)

data_pro<-read.table( 'climate.txt', sep="\t", header=T)

data_inpput<-xts(data_pro[, -1], order.by=as.Date(data_pro$Date))
other_para<- read.table( 'other_para_PET.txt', sep="\t", header=T)

data("constants")

constants$Elev<-other_para$Elev


constants$lat_rad<-other_para$lat_rad
constants$z<-other_para$z

constants$as<-0.12  # this value cite from Liu changming Acta Geographica sinica 2011 
constants$bs<-0.46

data_pet<-list(Date.daily=as.Date(data_pro$Date),J= yday(as.Date(data_pro$Date) ), 
               uz=data_inpput$uz, n =data_inpput$n, Tmin=data_inpput$Tmin, 
               Tmax =data_inpput$Tmax, RHmin=data_inpput$RHmin, RHmax=data_inpput$RHmax  )
PET_all <- ET.Penman(data_pet, constants, ts="daily",
                     solar="sunshine hours", wind="yes",
                     windfunction_ver = "1956", alpha = 0.23, z0 = 0.008)

#means1= as.xts(daily2annual(data_inpput$Pre, FUN=sum, na.rm=TRUE))

Pre<-as.matrix(data_inpput$Pre)
Runoff<-as.matrix(data_inpput$Runoff*24*3600/(1000*531))
Temp<-as.matrix(data_inpput$Ave_T)
PET<-as.matrix(PET_all$ET.Daily)
areas<-531 #the area of basin 

data(example_TUWmodel)


#-----------------------------
# Calibration module
#-----------------------------
#the range of parametrs 
paras_all<-list (  
  csf=c( 0.9,1.5), 
  ddf=c(0.0,5.0  ),
  tr=c(1.0,3.0  ),
  ts=c(-3.0,1.0 ), 
  tm=c(-2.0,2.0  ),
  lprat=c(0.0,1.0 ),
  fc=c(0,600),
  beta=c(0.0,20.0 ),
  k0=c(0.0,2.0  ),
  k1=c( 2.0,30.0 ),
  k2=c( 30.0,250.0  ),
  lsuz=c(1.0,100.0  ),
  cperc=c( 0.0,8.0  ),
  bmax=c(  0.0,30.0  ),
  croute=c(  0.0,50.0 )  
)

nn<-8400
Model.main<-function(paras){
  Qobs<-Runoff[1:nn]
  simDist1 <- HBV_E(prec=Pre, airt=Temp, ep=PET, area=areas/sum(areas),
                    param= paras,incon=c(60,0,1,2.5),iLength=nn)
  Qs<-simDist1$Q
  return(Qs)
}


control <- list(nseq=20,ndraw=50)
calbriate<-dreamCalibrate(
  FUN=Model.main,
  pars=paras_all,
  obs=Qobs,
  control=control
)

param_cal_out<- capture.output(summary(calbriate))


cat("param_cal_out", param_cal_out, file="param_cal_out.txt", sep="n", append=TRUE)
print(coef(calbriate)) 
#plot(calbriate)

#dd$lik.fun(coef(dd))
#
para_dream<-param_cal_out[52:66]
par_median_all<-par_2.5q_all<-par_97.5q_all<-NULL
for (ii in 1:15) {  
  pat <- "[-+.e0-9]*\\d"
  ac<-para_dream[[ii]]
  
  if ( ii <=11 & ii>=9) {
    par_median<-sapply(ac, function(x) strapply(x, pat, as.numeric)[[1]])[4]
    par_2.5q<-sapply(ac, function(x) strapply(x, pat, as.numeric)[[1]])[2]
    par_97.5q<-sapply(ac, function(x) strapply(x, pat, as.numeric)[[1]])[6]
  } else {
    par_median<-sapply(ac, function(x) strapply(x, pat, as.numeric)[[1]])[3]
    par_2.5q<-sapply(ac, function(x) strapply(x, pat, as.numeric)[[1]])[1]
    par_97.5q<-sapply(ac, function(x) strapply(x, pat, as.numeric)[[1]])[5]
    
  }
  
  par_median_all<-cbind(par_median_all,par_median)
  par_2.5q_all<-cbind(par_2.5q_all,par_2.5q)
  par_97.5q_all<-cbind(par_97.5q_all,par_97.5q)
  
}


sim_output_medain <- HBV_E(prec=Pre, airt=Temp, ep=PET, area=1,
                           param=as.numeric(par_median_all),incon=c(60,0,1,2.5),iLength=NULL)
sim_output_2.5q <- HBV_E(prec=Pre, airt=Temp, ep=PET, area=1,
                         param=as.numeric(par_2.5q_all),incon=c(60,0,1,2.5),iLength=NULL)
sim_output_97.5q <- HBV_E(prec=Pre, airt=Temp, ep=PET, area=1,
                          param=as.numeric(par_97.5q_all),incon=c(60,0,1,2.5),iLength=NULL)
sim_output<-cbind( runoff=as.numeric(sim_output_medain$Q), soil_moisture= as.numeric(sim_output_medain$moist), ETA=as.numeric(sim_output_medain$eta) )
rownames(sim_output)<- rownames(Runoff) 

write.table(sim_output, file = 'sim_output.txt', sep = ",", col.names = NA )



plot(as.Date(rownames(Runoff)), Runoff, type="l", ylab="Runoff [mm/day]", xlab='Year')
lines(as.Date(rownames(Runoff)), sim_output_medain$Q, col=2)

#lines(as.Date(rownames(Runoff)), sim_output_2.5q$Q, col = 'grey')
#lines(as.Date(rownames(Runoff)), sim_output_97.5q$Q, col = 'grey')
#polygon(c(as.Date(rownames(Runoff)), rev(as.Date(rownames(Runoff)))), c( sim_output_97.5q$Q, rev(sim_output_2.5q$Q)),
#     col = "grey30", border = NA)
legend("topleft", legend=c("Observations","Simulations"), col=c(1,2), lty=1, bty="n")



nse_value<-NSE(as.numeric(Runoff), as.numeric(sim_output_medain$Q) )

rmse_value<-rmse( as.numeric(sim_output_medain$Q),as.numeric(Runoff) )

r_value<-cor(as.numeric(Runoff), as.numeric(sim_output_medain$Q) )



#-----------------------------
# HBV-E main module
#-----------------------------


HBV_E <- function (prec, airt, ep, area=1, param=c(1.2,1.2,2,-2,0,0.9,100,3.3,0.5,9,105,50,2,10,26.5), incon=c(60,0,1,2.5), iLength=100) {
  nzones <- ifelse(is.vector(prec), 1, dim(prec)[2])
  iLength <- ifelse(is.null(iLength), length(prec)/nzones, iLength)
  if (nzones == 1) {
    
    parametri <- t(as.data.frame(t(param)))
    
    inconditions <- t(as.data.frame(t(incon)))
  } else if (nzones < 1) {
    cat("\nFormatting harddisk. Please smile!\n")
  } else {
    if (is.matrix(param)) {
      parametri <- param
    } else if (is.vector(param)) {
      parametri <- matrix(rep(param, nzones), ncol=nzones)
    }
    if (is.matrix(incon)) {
      inconditions <- incon
    } else if (is.vector(incon)) {
      inconditions <- matrix(rep(incon, nzones), ncol=nzones)
    }
  }
  storage.mode(prec) <- "double"
  storage.mode(airt) <- "double"
  storage.mode(ep) <- "double"
  storage.mode(area) <- "double"
  storage.mode(parametri) <- "double"
  storage.mode(inconditions) <- "double"
  output <- array(-777, dim=c(nzones, 20, iLength))
  storage.mode(output) <- "double"
  
  #dummy <- .Fortran("hbvmodel", iLength=as.integer(iLength), nzones=as.integer(nzones), area=area, param=parametri, incon=inconditions,
  #     prec=prec, airt=airt, ep=ep, output=output, PACKAGE="TUWmodel") 
  
  dummy <-  hbvmodel(iLength=as.integer(iLength), nzones=as.integer(nzones), area=area, param=parametri, incon=inconditions,
                     prec=prec, airt=airt, ep=ep, output=output  ) 
  
  #iLength=as.integer(iLength); nzones=as.integer(nzones); area=area; param=parametri; incon=inconditions;
  #prec=prec; airt=airt; ep=ep; output=output;
  
  
  names(dummy$param) <- c("SCF","DDF","Tr","Ts","Tm","LPrat","FC","BETA","k0","k1","k2","lsuz","cperc","bmax","croute")
  names(dummy$incon) <- c("SSM0","SWE0","SUZ0","SLZ0")
  dummy$qzones <- t(dummy$output[,1,])
  if (nzones > 1) {
    dummy$Q <- apply(dummy$qzones,1,weighted.mean,w=area)
  } else {
    dummy$Q <- dummy$qzones
  }
  dummy$swe <- t(dummy$output[,2,])
  dummy$melt <- t(dummy$output[,6,])
  dummy$q0 <- t(dummy$output[,7,])
  dummy$q1 <- t(dummy$output[,8,])
  dummy$q2 <- t(dummy$output[,9,])
  dummy$moist <- t(dummy$output[,3,])
  dummy$rain <- t(dummy$output[,4,])
  dummy$snow <- t(dummy$output[,5,])
  dummy$eta <- t(dummy$output[,10,])
  dummy$suz <- t(dummy$output[,11,])
  dummy$slz <- t(dummy$output[,12,])
  
  
  
  return(dummy)
}


hbvmodel<- function(iLength,nzones,area,param,incon, prec,airt,ep,output) { 
  
  
  imodincon=4; maxout=20; maxparam=15
  
  snowd<-matrix(nrow=iLength, ncol=nzones)
  
  #other temporal paras 
  #csf=ddf=tr=ts=meltt=lprat=LP=FC=beta=k0=k1=k2=lsuz=cperc=bmax=croute=sweq=NULL # sweq with a length of iLength 
  sweq=NULL 
  #sweq=rep(NULL, iLength)
  #not to define those variabes or parameter, just list how many parameters inovlved, would be helpful for users to understand. 
  #temp=precip=swe=NULL
  rain=snow=melt=NULL
  #moist=sd=
  dmoist=dq=eta=NULL
  #suz=slz=dquh=qg=q0=q1=q2=age==NULL
  qg=q0=q1=q2=age=NULL
  bql=NULL
  
  dquh=Q=assimilation=NULL
  
  # output[,,]=0
  # param<-t(param)
  
  #incon<-t(incon)
  for (izone in 1:nzones) {
    
    dmoist=dq=eta=NULL
    dmoist=dq=eta=NULL
    qg=q0=q1=q2=age=NULL
    bql=NULL
    dquh=Q=assimilation=NULL
    
    
    
    csf <- param[1,izone]
    ddf <- param[2,izone]
    tr <- param[3,izone]
    ts <- param[4,izone]
    meltt <- param[5,izone]
    
    lprat <- param[6,izone]
    FC <- param[7,izone]
    beta <- param[8,izone]
    
    LP <- lprat*FC
    
    k0 <- param[9,izone]
    k1 <- param[10,izone]
    k2 <- param[11,izone]
    lsuz <- param[12,izone]
    
    cperc <- param[13,izone]
    bmax <- param[14,izone]
    croute <- param[15,izone]
    
    ddfmin <- 1
    ddfmax <- 6
    ddfage <- ddf
    age<-0
    moist<-incon[1,izone]
    swe<-incon[2,izone]   
    suz<-incon[3,izone]   
    slz<-incon[4,izone] 
    sdold<-0 #
    
    if(area[izone]> 0) {
      for (i in 1:iLength){
        Q[i] <- 0
        dquh[i] <- 0
      }  
      aa <- area[izone]
      
      for (it in 1:iLength) {
        
        
        precip <- prec [it,izone]
        temp <- airt[it,izone]
        etp <- ep[it,izone]
        sd <- snowd[it,izone]
        
        
        if (temp<  -0.1 ) {
          etp <- 0
        }
        
        
        if (precip < -998.00) {
          for (ii in 1:12) {
            
            output[izone,ii,it] <- -999.99
          }
          
          
        } else {
          snowmod_out<-snowmod(csf,ddf,tr,ts,meltt,temp,precip,swe,  rain,snow,melt)
          
          csf=snowmod_out[,1];ddf=snowmod_out[,2];tr=snowmod_out[,3];ts=snowmod_out[,4];meltt=snowmod_out[,5];temp=snowmod_out[,6];precip=snowmod_out[,7];swe=snowmod_out[,8]; 
          rain=snowmod_out[,9];snow=snowmod_out[,10];melt=snowmod_out[,11];
          
          
          soilmoisture_out<- soilmoisture(rain,melt,etp,LP,FC,beta,dmoist,moist,  dq,eta)
          
          rain=soilmoisture_out[,1];melt=soilmoisture_out[,2];etp=soilmoisture_out[,3];LP=soilmoisture_out[,4];FC=soilmoisture_out[,5];beta=soilmoisture_out[,6];
          dmoist=soilmoisture_out[,7];moist=soilmoisture_out[,8];  dq=soilmoisture_out[,9];eta=soilmoisture_out[,10]
          
          
          respfunc_out<- respfunc(dq,k0,lsuz,k1,k2,cperc,bmax,croute,suz,slz,bql,dquh,qg,q0,q1,q2)
          
          dq=respfunc_out[,1];k0=respfunc_out[,2];lsuz=respfunc_out[,3];k1=respfunc_out[,4];k2=respfunc_out[,5];cperc=respfunc_out[,6];bmax=respfunc_out[,7];croute=respfunc_out[,8];
          suz=respfunc_out[,9];slz=respfunc_out[,10];bql=respfunc_out[,11];dquh=respfunc_out[,12];qg=respfunc_out[,13];q0=respfunc_out[,14];q1=respfunc_out[,15];q2=respfunc_out[,16]
          sum_ml=respfunc_out[,17]
          
          
          
          
          for ( irf in 1:bql) {
            if ((it+irf-1)<=iLength) {      
              Q[it+irf-1]=Q[it+irf-1]+dquh[irf] } #DAMP LAG route in rivers 
          }
          
          
          output[izone,1,it]=Q[it] 
          output[izone,2,it]=swe
          sweq[it]=swe     
          output[izone,3,it]=moist
          output[izone,4,it]=rain
          output[izone,5,it]=snow
          output[izone,6,it]=melt
          output[izone,7,it]=q0[1]
          output[izone,8,it]=q1[1]
          output[izone,9,it]=q2[1]
          output[izone,10,it]=eta
          output[izone,11,it]=suz[1]    
          output[izone,12,it]=slz [1]   
          
          
        }
      }
      
      
    } else  {
      for (it in 1:iLength) {
        for (ii in 1:12) {
          
          output[izone,ii,it] <- 0
        }
        
      }
      
    }
    
    
  }
  
  results<-list(iLength=iLength,nzones=nzones,area=area,param=param,incon=incon, prec=prec,airt=airt,ep=ep,output=output)
  return(results)
}



respfunc <- function (dq,k0,lsuz,k1,k2,cperc,bmax,croute,suz,slz,bql,dquh,qg,q0,q1,q2) {
  
  maxday<-15000
  dt=1.0
  rat=1.0
  suzold<-suz+rat*dq
  slzold<-slz+(1-rat)*dq
  slzin<-cperc
  
  if(suzold <0) suzold=0
  if(slzold < 0) slzold=0
  
  #     --- 1st storage ---
  if (suzold>lsuz){
    q0=(suzold-lsuz)/k0*exp(-1/k0)
    if (q0<0){
      q0 <- 0
    }
    if (q0>(suzold-lsuz)){
      q0 <- suzold-lsuz
    }
  } else {
    q0 <- 0
  }
  
  suzold <- suzold-q0
  q1 <- -slzin+(slzin+suzold/k1)*exp(-1/k1)
  
  if (q1<0 ){
    q1 <- 0
  }
  suz <- suzold-q1-slzin
  if (suz<0){
    suz <- 0
    slzin <- suzold
  }
  
  #     --- 2rd storage ---
  q2=slzin-(slzin-slzold/k2)*exp(-1/k2)
  if (q2<0)     q2 <- 0  
  
  slz <- slzold-q2+slzin
  if (slz<0){
    slz <- 0
    q2 <- slzold+slzin
  }
  
  
  qg <- q0+q1+q2
  
  #     --- transformation function ---
  if ((bmax-croute*qg)>1){
    bq <- bmax-croute*qg
    bql <-  round(bq,0)
    sum_ml <- 0
    for ( j in 1:bql){
      
      if (j<=bql/2){
        dquh[j]=((j-0.5)*4*qg)/(bql*bql*1)
        
      } else if (abs(j-(bql/2+0.5))<0.1) {
        dquh[j] <-((j-0.75) *4*qg)/(bql*bql*1)
      } else {
        dquh[j] <-((bql-j+0.5)*4*qg)/(bql*bql*1)
      }
      sum_ml <- sum_ml+dquh[j]
    }
    
    
  } else {
    bql <- 1
    dquh[1] <- qg
    sum_ml <- qg
  }
  
  respfunc_out<-cbind(dq,k0,lsuz,k1,k2,cperc,bmax,croute,suz,slz,bql,dquh,qg,q0,q1,q2,sum_ml)
  return (respfunc_out) 
  
  
}



soilmoisture <- function(rain,melt,etp,LP,FC,beta,dmoist,moist, dq,eta) { 
  
  xx=NULL
  #     --- soil mositure accounting ----
  moistold <- moist
  dq <-((moistold/FC)**beta)*(rain+melt)     
  if (dq>(rain+melt))  dq <- rain+melt
  dmoist <- rain+melt-dq
  if (dmoist<0)   dmoist <- 0
  
  
  moist <- moistold+dmoist
  if (moist>FC){
    dq <- (moist-FC)+dq
    moist <- FC
  }
  #     --- calculate evapotranspiration ---
  if (moist<LP){
    eta <- moist*etp/LP
    if (eta>etp){
      eta <- etp
    }
  } else {
    eta <- etp
  }
  if (eta<0)   eta <- 0   
  
  #     --- substract eta of soilmoisture ---
  xx <- moist
  moist <- moist-eta
  if (moist<0){
    eta <- xx
    moist <- 0
  }
  
  soilmoisture_out<-cbind(rain,melt,etp,LP,FC,beta,dmoist,moist, dq,eta) 
  return(soilmoisture_out)
  
}



snowmod<- function (csf,ddf,tr,ts,melttemp,temp,precip,swe,rain,snow,melt) { 
  
  if (temp<ts){
    snow <- precip
  } else if (temp>tr) {
    snow <- 0
  } else {
    snow <- precip*abs(temp-tr)/abs(tr-ts)
  }
  
  rain <- precip-snow
  
  
  melt <- (temp-melttemp)*ddf
  if (melt < 0){
    melt <- 0
  }
  
  sweold <- swe
  swe <- sweold+csf*snow-melt
  if (swe< 0.0001){
    swe <- 0
    melt <- sweold+csf*snow
    if (melt<0)    melt <- 0  
  }
  
  
  snowmod_out<-cbind(csf,ddf,tr,ts,melttemp,temp,precip,swe,rain,snow,melt )
  
  return(snowmod_out)
  
}



