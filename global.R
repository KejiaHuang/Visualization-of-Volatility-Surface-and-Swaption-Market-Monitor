
library(shiny)
library(plotly)
library(lubridate)
library(shinydashboard)
setwd("D:\\FE\\800\\666shiny\\version_2.0")
Vols_data<-read.csv("Swaption_vols.csv")
Strikes_data<-read.csv("Swaption_strikes.csv")

# define matrix dimension
row_num <- dim(Vols_data)[1]
col_num <- dim(Vols_data)[2]

# define function to get upper limit
up_limit <- function(vol){
  # test data
  # vol <- Vols_data[,1]
  mean <- mean(as.numeric(vol))
  sd  <-sd(as.numeric(vol))
  up_result <- mean + 1.96 * sd 
  return(up_result)
}


# define function to get lower limit
low_limit <- function(vol){
  # test data
  # vol <- Vols_data[,1]
  mean <- mean(as.numeric(vol))
  sd  <-sd(as.numeric(vol))
  low_result <- mean - 1.96 * sd 
  return(low_result)
}

# use apply function to build up_matrix and low_matrix

cut_matrix  <- c(1:row_num)
cut_matrix = cut_matrix[cut_matrix %% 30 == 0]
up_data_sum <- matrix(nrow = 8, ncol = 271)
low_data_sum <- matrix(nrow = 8, ncol = 271)
j=1
for(i in cut_matrix){
  up_data  <- apply( Vols_data[(i-29):(i),-1], 2, up_limit  )
  low_data <- apply( Vols_data[(i-29):(i),-1], 2, low_limit )  
  up_data_sum[j,]  <- matrix(c(as.character(Vols_data[i,1]),up_data ), nrow =1)
  low_data_sum[j,] <- matrix(c(as.character(Vols_data[i,1]),low_data), nrow =1)
  j = j + 1
  
}




##################################################################

library(plotly)

load_data<-function(date,Vols_data,Strikes_data){
  # vol_data[vol_data$Date=="2016/1/1",]
  # date<-"2016/1/1"; Vols_data<-vol_data; Strikes_data<-K_data
  vols<-Vols_data[as.Date(Vols_data$Date)==as.Date(date),]
  vols<-vols[,-1]
  vols<-matrix(unlist(vols),nrow = 15, ncol = 18, byrow = FALSE )#
  vols<-t(vols)
  strikes<-Strikes_data[as.Date(Strikes_data$Date)==as.Date(date),]
  strikes<-strikes[,-1]
  strikes<-matrix(unlist(strikes),nrow = 15, ncol = 18, byrow = FALSE )#
  strikes<-t(strikes)
  return(list(vols,strikes))
  
}



SABRVol<-function(alpha,beta,rho,nu,T,f,K){
  
  if(f==K){
    
    Denom<- f^(1-beta)
    Term1<- 1+(((1-beta)^2)/24*(alpha^2)/(f^(2-2*beta))+1/4*(rho*beta*nu*alpha)/(f^(1-beta))+(2-3*rho^2)/24*nu^2)*T
    return(alpha/Denom*Term1)
    
  }
  
  else{
    z<-nu/alpha*((f*K)^(1-beta)/2)*log(f/K)
    X<-function(z){
      return(log((sqrt(1-2*rho*z+z^2)+z-rho)/(1-rho)))
    }
    Denom <-(f*K)^((1-beta)/2)*(1+((1-beta)^2)/24*log(f/K)^2+((1-beta)^4)/1920*log(f/K)^4)
    Term1<- 1+(((1-beta)^2)/24*alpha^2/((f*K)^(1-beta))+1/4*(rho*beta*nu*alpha)/((f*K)^((1-beta)/2))+(2-3*rho^2)/24*nu^2)*T
    return(alpha/Denom*Term1*z/X(z))
  }
  
}
FitSABR<-function(beta,fATM,K,t,obvols,
                  init.values = c(0.0, 0, 0),
                  lower.bound = c(0.00, -1, 0.0001),
                  upper.bound = c(1, 1, 5)){
  #beta=0.5;  fATM=K_data[1,]; K=K_data[1,]; t=Ex[,1]; obvols=test_data;
  obj<-function(parm,beta,fATM,K,t,obvols){
    #alpha<-0.11;rho=-0.2;nu=2;
    alpha<-parm[1]
    rho<-parm[2]
    nu<-parm[3]
    ESvols<-SABRVol(alpha=alpha,beta=beta,rho=rho,nu=nu,T=t,f=fATM,K=K)
    return(sum((obvols-ESvols)^2))
  }
  #obj(c(0.11,-0.2,2),beta,fATM,K,t,obvols)
  #
  opt<-nlminb(start=init.values,
              objective = obj,
              lower = lower.bound,
              upper = upper.bound,
              beta=beta,
              fATM=fATM,
              K=K,
              t=t,
              obvols=obvols)
  parms<-opt$par
  names(parms) <- c("alpha", "rho", "nu")
  obj <- opt$objective
  # opt<-optim(par=c(0, 0, 1),obj,beta=beta,fATM=fATM,K=K,t=t,obvols=obvols)
  
  ESvols<-SABRVol(alpha=parms[1],beta=0.5,rho=parms[2],nu=parms[3],T=t,f=fATM,K=K)
  #ESvols-obvols
  return(list(parms=parms,  estvols=ESvols,  obsvols=obvols,
              obj=obj))
}

tenor<- c(1,2,3,4,5,6,7,8,9,10,12,15,20,25,30)
Ex <-c(1/12,3/12,6/12,9/12,1,2,3,4,5,6,7,8,9,10,15,20,25,30)
expiry<-c(1/12,3/12,6/12,9/12,1,2,3,4,5,6,7,8,9,10,15,20,25,30)

sabr_row <- function(x){
  vols<-matrix(unlist(Vols_data[x,-1])/100,nrow = 15, ncol = 18 , byrow = FALSE)
  strikes<-matrix(unlist(Strikes_data[x,-1])/100,nrow = 15, ncol = 18, byrow = FALSE )
  SABR_Parms<-matrix(nrow=nrow(vols),ncol=3)
  Par <- unlist(lapply( 1:15, function(x) 
    FitSABR(0.5,strikes[x,],strikes[x,],Ex,vols[x,])$parms))
  Parms <- matrix(c(as.character(Vols_data[x,1]),as.numeric(format(round(Par,4)))),nrow=1, byrow = TRUE)
  return(Parms)
}

row_index <- c(1:row_num)
row_index = row_index[row_index %% 30 == 0]
final <- lapply(row_index,sabr_row)
Parms <- matrix(unlist(final),nrow = 8, byrow = TRUE)
tenor_plus<-c()
for(i in tenor){
  # i = tenor[1]
  tenor_plus<-c(tenor_plus,rep(i,18))
}
