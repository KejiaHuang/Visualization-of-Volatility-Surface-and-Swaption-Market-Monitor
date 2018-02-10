##### global--------
options(warn = -1) 
Sys.setlocale("LC_TIME", "English")
library(shiny)
library(plotly)
library(shinydashboard)
library(lubridate)

setwd("D:\\FE\\800\\666shiny\\version_final_modified")

Vols_data<-read.csv("Swaption_vols.csv")
Strikes_data<-read.csv("Swaption_strikes.csv")
yield_data<-read.csv("yield_curve.csv")

tenor<- c(1,2,3,4,5,6,7,8,9,10,12,15,20,25,30)
Ex <-c(1/12,3/12,6/12,9/12,1,2,3,4,5,6,7,8,9,10,15,20,25,30)
expiry<-c(1/12,3/12,6/12,9/12,1,2,3,4,5,6,7,8,9,10,15,20,25,30)


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
#---Q1
FitSABR<-function(beta,fATM,K,t,obvols,
                  init.values = c(10, 0, 3),
                  lower.bound = c(0.00, -1, 0.0001),
                  upper.bound = c(100, 1, 100)){
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
  # names(parms) <- c("alpha", "rho", "nu")
  obj <- opt$objective
  # opt<-optim(par=c(0, 0, 1),obj,beta=beta,fATM=fATM,K=K,t=t,obvols=obvols)
  
  ESvols<-SABRVol(alpha=parms[1],beta=0.5,rho=parms[2],nu=parms[3],T=t,f=fATM,K=K)
  #ESvols-obvols
  return(list(parms=parms,  estvols=ESvols,  obsvols=obvols,
              obj=obj))
}


SABRVol<-function(alpha,beta,rho,nu,T,f,K){
  
  # if(f==K){
    
    Denom<- f^(1-beta)
    Term1<- 1+(((1-beta)^2)/24*(alpha^2)/(f^(2-2*beta))+1/4*(rho*beta*nu*alpha)/(f^(1-beta))+(2-3*rho^2)/24*nu^2)*T
    return(alpha/Denom*Term1)
    
  # }
  
  # else{
  #   z<-nu/alpha*((f*K)^(1-beta)/2)*log(f/K)
  #   X<-function(z){
  #     return(log((sqrt(1-2*rho*z+z^2)+z-rho)/(1-rho)))
  #   }
  #   Denom <-(f*K)^((1-beta)/2)*(1+((1-beta)^2)/24*log(f/K)^2+((1-beta)^4)/1920*log(f/K)^4)
  #   Term1<- 1+(((1-beta)^2)/24*alpha^2/((f*K)^(1-beta))+1/4*(rho*beta*nu*alpha)/((f*K)^((1-beta)/2))+(2-3*rho^2)/24*nu^2)*T
  #   return(alpha/Denom*Term1*z/X(z))
  # }
  
}
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
# monitor new
load_data30<-function(date,Vols_data,Strikes_data){
  # vol_data[vol_data$Date=="2016/1/1",]
  # date<-"2016/1/1"; Vols_data<-vol_data; Strikes_data<-K_data
  vols<-Vols_data[as.Date(Vols_data$Date)==as.Date(date),]
  location <- which(as.Date(Vols_data$Date)==as.Date(date))
  vols_sum<-Vols_data[(location-30):(location-1),-1]/100
  strikes_sum<- Strikes_data[(location-30):(location-1),-1]/100
  
  return(list(vols_sum,strikes_sum))
  
}
row_Sabr<-function(i,vols_data30,strikes_data30){
  row_data_vol<-vols_data30[i,]
  row_data_K<-strikes_data30[i,]
  vols<-matrix(unlist(row_data_vol),nrow = 15, ncol = 18, byrow = FALSE )#
  strikes<-matrix(unlist(row_data_K),nrow = 15, ncol = 18, byrow = FALSE )#
  parms<-matrix(nrow=15,ncol=3)
  SABR_vols<-t(apply(matrix(c(1:15),nrow=15),1,apply_row_Sabr<-function(i){
    parms[i,]<-FitSABR(0.5,strikes[i,],strikes[i,],Ex,vols[i,])$par
    SABR_vols<-matrix(nrow=15,ncol=18)
    SABR_vols[i,]<-SABRVol(parms[i,1],0.5,parms[i,2],parms[i,3],Ex,strikes[i,],strikes[i,])
    return(SABR_vols[i,])
  }))
  SABR_vols<-split(t(SABR_vols), rep(1))[[1]]
  return(SABR_vols)
}

up_limit_new <- function(vol,c){
  # test data
  # vol <- Vols_data[,1]
  alpha<-(1-c)/2
  ci<-abs(qnorm(alpha))
  mean <- mean(as.numeric(vol))
  sd  <-sd(as.numeric(vol))
  up_result <- mean + ci * sd 
  return(up_result)
}

low_limit_new <- function(vol,c){
  # test data
  # vol <- Vols_data[,1]
  # c<-0.95
  alpha<-(1-c)/2
  ci<-abs(qnorm(alpha))
  mean <- mean(as.numeric(vol))
  sd  <-sd(as.numeric(vol))
  low_result <- mean - ci * sd 
  return(low_result)
}

yield_curve<-function(date,yield_data){
  if(weekdays(date)=="Saturday"){
    date<-as.Date(date)+2
  }
  else if(weekdays(date)=="Sunday"){
    date<-as.Date(date)+1
  }
  else{}
  yield_c<-yield_data[as.Date(yield_data$Date)==as.Date(date),]
  yield_c<-yield_c[,-1]
  return(yield_c)
}


BS_Swaption<-function(S,K,ts,T,sigma,df_data){
  #S<-0.12;K<-0.12;ts<-1;T<-1*12;sigma<-0.5,
  dt<-0.5
  time_sires<-seq(from = 6 , to = 360 ,by = 6)
  T_mo<-T*12
  # time_sires<-c(time_sires)
  # loc_ts<-which(time_sires==ts)
  loc_T<-which(time_sires==T_mo)
  df<-df_data[1:loc_T]
  DV01<-sum(df)*dt
  d1<-(log(S/K)+sigma^2*ts/2)/(sigma*sqrt(ts))
  d2<-d1-sigma*sqrt(ts)
  swaption<-DV01*(S*pnorm(d1)-K*pnorm(d2))
  return(swaption)
}

calculate_bs <- function(vol,strike,date,yield_data){
  #vol<-today_upper_vols; strike<-today_upper_k
  price<-matrix(nrow=15,ncol=18)
  tenor<- c(1,2,3,4,5,6,7,8,9,10,12,15,20,25,30)
  Ex <-c(1/12,3/12,6/12,9/12,1,2,3,4,5,6,7,8,9,10,15,20,25,30)
  yield_curve_temp <- yield_curve(date,yield_data)
  for(i in 1:15){
    for(j in 1:18){
      if(!is.na(vol[i,j])){
        price[i,j]<-BS_Swaption(strike[i,j],strike[i,j],Ex[j],tenor[i],vol[i,j],yield_curve_temp)
      }
    }
  }
  return(price)
}

##### ui    -------     -----
## app.R ##
library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(title = "volitility surface"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("SABR model", tabName = "dashboard", icon = icon("area-chart")),
      menuItem("Calibration", tabName = "widgets", icon = icon("dashboard")),
      #menuItem("Parameters", tabName = "parms", icon = icon("calendar")),
      #menuItem("SABR calculator",tabName="P2",icon=icon("dashboard")),
      #menuItem("Monitor", tabName = "moditor", icon = icon("binoculars")),
      menuItem("SABR Monitor",tabName="monitor_new" ,icon=icon("binoculars")),
      menuItem("Info", tabName = "info", icon = icon("wechat"))
    )
  ),
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "dashboard",
              fluidRow(
                
                box(
                  title = "Calendar",
                  dateInput("date1", "Choose a date to display:", 
                            value = "2016-4-20", min = "2016-1-1", max = "2017-1-1"),
                  width = 3,
                  solidHeader = TRUE,
                  status = "primary"
                ),
                box(
                    title = "Inputs",
                      checkboxInput("flag1", "Observed Vol", TRUE),
                     checkboxInput("flag2", "SABR Vol", TRUE),
                     solidHeader = TRUE,
                     status = "warning",
                    width=3
                     ),
                box( plotlyOutput("plot1" ),width=9)
                
               
              )
      ),
      
      
   
      
      # Second tab content
      tabItem(tabName = "widgets",
              fluidRow(
                box(
                  title = "Calendar",
                  dateInput("date2", "Choose a date to display:", 
                            value = "2016-4-20", min = "2016-1-1", max = "2017-1-1"),
                  width = 3,
                  solidHeader = TRUE,
                  status = "primary"
                ),
                box( plotOutput("plot2" ),width = 12,hight = 15),
                box(
                  title="Inputs",
                  width = 5,selectInput("select", label = h3("Tenor"),
                                        choices = list("1Y" = 1, "2Y" = 2,"3Y" = 3,"4Y" = 4,"5Y" = 5,
                                                       "6Y" = 6, "7Y" = 7,"8Y" = 8,"9Y" = 9,"10Y" = 10,
                                                       "12Y" = 11, "15Y" = 12,"20Y" = 13,"25Y" = 14,"30Y" = 15,value=12
                                        )),
                  solidHeader = TRUE,
                  status = "warning")
                ,
                valueBox(0.5, h2("Beta"),width = 3, icon = icon("cogs")),
                valueBoxOutput("Alpha",width = 3),
                valueBoxOutput("Rho",width = 3),
                valueBoxOutput("Nu",width = 3),
                box( plotlyOutput("plot4" ),width = 12)
                
              )
      ),
      # ##### parameter 
      # tabItem(tabName = "parms",
      #         fluidRow(
      # 
      #           box(
      #             title = "Calendar",
      #             dateInput("date4", "Choose a date to display:",
      #                       value = "2016-4-20", min = "2016-1-1", max = "2017-1-1"),
      #             width = 3,
      #             solidHeader = TRUE,
      #             status ="primary"
      #           ),
      #           #1,2,3,4,5,6,7,8,9,10,12,15,20,25,30
      #           box(
      #               title="Inputs",
      #               width = 3,selectInput("select", label = h3("Tenor"),
      #                            choices = list("1Y" = 1, "2Y" = 2,"3Y" = 3,"4Y" = 4,"5Y" = 5,
      #                                           "6Y" = 6, "7Y" = 7,"8Y" = 8,"9Y" = 9,"10Y" = 10,
      #                                           "12Y" = 11, "15Y" = 12,"20Y" = 13,"25Y" = 14,"30Y" = 15,value=12
      #                                           )),
      #                solidHeader = TRUE,
      #                status = "warning")
      #           ,
      #           valueBox(0.5, h2("Beta"),width = 3, icon = icon("cogs")),
      #           valueBoxOutput("Alpha",width = 3),
      #           valueBoxOutput("Rho",width = 3),
      #           valueBoxOutput("Nu",width = 3),
      #           box( plotlyOutput("plot4" ),width = 9))
      # 
      # 
      # 
      #         
      # ),
     tabItem(tabName="P2",
              fluidRow(
                
                box(
                  title = "Calendar",
                  dateInput("datex", "Choose a date to display:",
                            value = "2016-4-21", min = "2016-1-1", max = "2017-1-1"),
                  width = 3,
                  solidHeader = TRUE,
                  status ="primary"
                ),
                #1,2,3,4,5,6,7,8,9,10,12,15,20,25,30
                box(
                  title="Inputs",
                  width = 3,selectInput("selectx", label = h3("Tenor"),
                                        choices = list("1Y" = 1, "2Y" = 2,"3Y" = 3,"4Y" = 4,"5Y" = 5,
                                                       "6Y" = 6, "7Y" = 7,"8Y" = 8,"9Y" = 9,"10Y" = 10,
                                                       "12Y" = 11, "15Y" = 12,"20Y" = 13,"25Y" = 14,"30Y" = 15
                                        )),
                  solidHeader = TRUE,
                  status = "warning")
                ,
                box(
                sliderInput("betax",
                            "Beta",
                            min = 0.0001,
                            max = 1,
                            value = 0.5,
                            step=0.001),
                sliderInput("alphax",
                            "Alpha",
                            min = 0.0001,
                            max = 1,
                            value = 0.056,#parsx()[1],
                            step=0.001),
                sliderInput("rhox",
                            "Rho",
                            min = -1,
                            max = 1,
                            value = -1,#parsx()[2],
                            step = 0.01),
                sliderInput("nux",
                            "Nu",
                            min = 0.0001,
                            max = 5,
                            value = 0.1,#parsx()[3],
                            step=0.001)),
                box( plotlyOutput("plotx" ),width = 9))),
              
    
      # monitor
      tabItem(tabName = "moditor",
              fluidRow(
                
                box(
                  title = "Calendar",
                  dateInput("date3", "Choose a date to display:", 
                            value = "2016-4-21", min = "2016-1-1", max = "2017-1-1"),
                  width = 3,
                  solidHeader = TRUE,
                  status ="primary"
                ),
                box( plotlyOutput("plot3" ),width = 9),
                box(plotOutput("plotgg"),width=12)
                
              
                
                
              )
      ),
     #### monitor_new
     tabItem(tabName = "monitor_new",
             fluidRow( 
               box(
                    title = "Calendar",
                    dateInput("date_new", "Choose a date to display:", 
                    value = "2016-4-21", min = "2016-1-1", max = "2017-1-1"),
               width = 4,
               solidHeader = TRUE,
               status ="primary"
             ),
             box(title = "Notional",textInput("notional","Notional", value = "10000",width="800px", placeholder = NULL),
                 width=4,
                 solidHeader = TRUE,
                 status ="primary"),
             valueBoxOutput("total_notional",width = 4),
             # valueBoxOutput("long_amount",width = 3),
             # valueBoxOutput("short_amount",width = 3),
             # valueBoxOutput("total_amount",width = 3),
             
             box( plotlyOutput("plot_new" ),width = 12),
             # box(plotOutput("plotgg_new"),width=9),
             box(plotOutput("plotgg_price"),width=10,hight="120%"),
             box(sliderInput("c",
                             "Confidence interval",
                             min = 0.750,
                             max = 0.999,
                             value = 0.950,
                             step=0.001),width = 2),
             # valueBoxOutput("profit",width = 3),
             valueBoxOutput("rate",width = 4),
             valueBoxOutput("red",width = 4),
             valueBoxOutput("green",width = 4),
             valueBoxOutput("profit",width = 4),
             valueBoxOutput("r_profit",width = 4),
             valueBoxOutput("g_profit",width = 4),
             box(plotlyOutput("yield_c"),width=9),
             box(plotOutput("plotgg_return"),width=9),
             box(title = "Visualization tools",
                 sliderInput( "return1",
               "Lower limit of return(%)",
               min = 0,
               max = 2,
               value = 1,
               step=0.01),
               sliderInput(
                 "return2",
                 "Upper limit of return(%)",
                 min = 0,
                 max = 2,
                 value = 1,
                 step=0.01),
               width=3
               
             )
             
             )),
           
             
             
             
            
      
      
      #info
      tabItem(tabName = "info",
              fluidRow(
                infoBox("Kejia Huang","khuang10@stevens.edu",icon = icon("credit-card") ),
                
                infoBox("Libo Duan","lduan1@stevens.edu",icon = icon("credit-card") ),
                
                infoBox("Zenghui Liu",  "zliu45@stevens.edu",icon = icon("credit-card") )
              )
              )
    )
  )
)
##### server-------
server <- function(input, output) {
##### reactive ---
  volatiltiy <- reactive({
    #date<-"2016/4/23"
    date <- input$date1
    if (weekdays(as.Date(date))=="Sunday"){
      date = as.Date(date)-2
    }
    if (weekdays(as.Date(date))=="Saturday"){
      date = as.Date(date)-1
    }

    data<-load_data(date,Vols_data,Strikes_data)
    vol_data<-t(as.data.frame(data[1]))/100
    K_data<-t(as.data.frame(data[2]))/100
    
    SABR_Parms<-matrix(nrow=nrow(vol_data),ncol=3)
    SABR_Volatilities<-matrix(nrow=nrow(vol_data),ncol=ncol(vol_data))
    for(i in 1:15){
      SABR<-FitSABR(0.5,K_data[i,],K_data[i,],Ex,vol_data[i,])
      SABR_Parms[i,]<-SABR$par
      SABR_Volatilities[i,]<-SABR$estvols
      #plot(x=Ex[i,],y=vol_data[i,],type="l",xlab="expiry",ylab="vols",main=paste0("X year into ",tenor[i]," year swaption"))
      #lines(x=Ex[i,],y=SABR$estvols,type="l",col=2)
    }
    
    return(list(SABR_Volatilities*100,vol_data*100))
    
  })
  sub_vol <- reactive({
    #date<-"2016/4/21"
    date2 <- input$date2
    if (weekdays(as.Date(date2))=="Sunday"){
      date2 = as.Date(date2)-2
    }
    if (weekdays(as.Date(date2))=="Saturday"){
      date2 = as.Date(date2)-1
    }
    data<-load_data(date2,Vols_data,Strikes_data)
    vol_data<-t(as.data.frame(data[1]))/100
    K_data<-t(as.data.frame(data[2]))/100
    
   
    
    return(list(vol_data,K_data))
    
  })
  sub_vol_parms <- reactive({
    #date<-"2016/4/21"
    date4 <- input$date2
    if (weekdays(as.Date(date4))=="Sunday"){
      date4= as.Date(date4)-2
    }
    if (weekdays(as.Date(date4))=="Saturday"){
      date4 = as.Date(date4)-1
    }
    data<-load_data(date4,Vols_data,Strikes_data)
    vol_data<-t(as.data.frame(data[1]))/100
    K_data<-t(as.data.frame(data[2]))/100
    
    
    
    return(list(vol_data,K_data))
    
  })
  sub_vol_parmsx <- reactive({
    #date<-"2016/4/21"
    date4 <- input$datex
    if (weekdays(as.Date(date4))=="Sunday"){
      date4= as.Date(date4)-2
    }
    if (weekdays(as.Date(date4))=="Saturday"){
      date4 = as.Date(date4)-1
    }
    data<-load_data(date4,Vols_data,Strikes_data)
    vol_data<-t(as.data.frame(data[1]))/100
    K_data<-t(as.data.frame(data[2]))/100
    
    
    
    return(list(vol_data,K_data))
    
  })
  moni_vol <- reactive({

    location <-which(as.Date(Vols_data[,1]) == as.Date(input$date3))
    i0 <- floor(location / 30)
    Parms_matrix <- matrix(as.numeric(Parms[i0,-1]),nrow = 15, ncol = 3, byrow =  TRUE)
    K <- matrix(as.numeric(Strikes_data[location,-1]/100),nrow=15,ncol=18, byrow = FALSE)
    
    vol_matrix <- matrix(nrow=15,ncol=18)
    for(i in 1:nrow(vol_matrix)){
      vol_matrix[i,]<-SABRVol(Parms_matrix[i,1],beta=0.5,Parms_matrix[i,2],Parms_matrix[i,3],
                              Ex,K[i,],K[i,])
    }
    ob_vols<-matrix(as.numeric(Vols_data[location,-1]/100),nrow=15,ncol=18, byrow = FALSE)
    up_mark  <- vol_matrix > matrix(as.numeric(up_data_sum[i0,-1] )/100,nrow=15,ncol=18,byrow = FALSE)
    low_mark <- vol_matrix < matrix(as.numeric(low_data_sum[i0,-1])/100,nrow=15,ncol=18,byrow = FALSE)
    upper_matrix<-matrix(nrow=15,ncol=18)
    lower_matrix<-matrix(nrow=15,ncol=18)
    for(i in 1:nrow(vol_matrix)){
      for(j in 1:ncol(vol_matrix)){
        if(up_mark[i,j]){
          upper_matrix[i,j]=vol_matrix[i,j]
        }
        else{}
        if(low_mark[i,j]){
          lower_matrix[i,j]=vol_matrix[i,j]
        }
        else{}
      }
    }
    upper <- round(upper_matrix*100,digits = 2)
    lower <- round(lower_matrix*100,digits = 2)
    
    upper_limit <- matrix(as.numeric(up_data_sum[i0,-1] )/100,nrow=15,ncol=18,byrow = FALSE)
    lower_limit <- matrix(as.numeric(low_data_sum[i0,-1])/100,nrow=15,ncol=18,byrow = FALSE)
    upper_limit <- round(upper_limit*100,digits = 1)
    lower_limit <- round(lower_limit*100,digits = 1)
    EX_char <- c(1/12,3/12,6/12,9/12,1,2,3,4,5,6,7,8,9,10,15,20,25,30)
    EX_char <- as.factor(EX_char)
    Ex_char_plus <- rep(EX_char,15)
    #tenor_char <-c("1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y","12Y","15Y","20Y","25Y","30Y")
    tenor_char <-c(1,2,3,4,5,6,7,8,9,10,12,15,20,25,30)
    tenor_char <- as.factor(tenor_char)
    tenor_char_plus<-c()
    for(i in tenor_char){
      tenor_char_plus<-c(tenor_char_plus,rep(i,18))
    }
    upper_plot<-split(t(upper), rep(1))[[1]]
    lower_plot<-split(t(lower), rep(1))[[1]]
    upper_lower<-c(upper_plot,lower_plot)
    upper_lower<-as.data.frame(upper_lower)
    Ex_char_plus<-as.data.frame(rep(Ex_char_plus,2))
    
    tenor_char_plus<-as.data.frame(rep(tenor_char_plus,2))
    col_flag <- as.data.frame(c(rep(1,270),rep(2,270)))
    # levels(col_flag)[levels(col_flag)=="1"] <- "Upper"
    # levels(col_flag)[levels(col_flag)=="2"] <- "Lower"
    
    upper_lower_plot <- cbind(upper_lower,Ex_char_plus,tenor_char_plus,col_flag)
    upper_lower_plot <- as.data.frame(upper_lower_plot)
    colnames(upper_lower_plot) <- c("df_upperlower","df_Ex","df_tenor","df_color")
    
  return(list(ob_vols,upper_matrix,lower_matrix,upper_lower_plot))
  
  })
  moni_new <- reactive({
    
    #date<-"2016/7/14"
    date <- input$date_new
    c <- input$c
    data<-load_data(date,Vols_data,Strikes_data)
    vol_data_today<-t(as.data.frame(data[1]))/100
    K_data_today<-t(as.data.frame(data[2]))/100
    tenor<- c(1,2,3,4,5,6,7,8,9,10,12,15,20,25,30)
    Ex <-c(1/12,3/12,6/12,9/12,1,2,3,4,5,6,7,8,9,10,15,20,25,30)
    expiry<-c(1/12,3/12,6/12,9/12,1,2,3,4,5,6,7,8,9,10,15,20,25,30)
    
    data_30<-load_data30(date,Vols_data ,Strikes_data)
    vols_data30<-data_30[[1]]
    strikes_data30<-data_30[[2]]
    SABR_30<-t(apply(matrix(c(1:30),nrow=30),1,row_Sabr,vols_data30,strikes_data30))
    #input c
    #c<-0.999
    upper_limit <- apply(SABR_30,2,up_limit_new,c)
    lower_limit <- apply(SABR_30,2,low_limit_new,c)
    upper_limit<-matrix(unlist(upper_limit),nrow = 15, ncol = 18, byrow = FALSE )#
    lower_limit<-matrix(unlist(lower_limit),nrow = 15, ncol = 18, byrow = FALSE )#
    
    upper_marks<-vol_data_today>upper_limit
    lower_marks<-vol_data_today<lower_limit
    upper_matrix<-matrix(nrow=15,ncol=18)
    lower_matrix<-matrix(nrow=15,ncol=18)
    
    for(i in 1:nrow(vol_data_today)){
      for(j in 1:ncol(vol_data_today)){
        if(upper_marks[i,j]){
          upper_matrix[i,j]=vol_data_today[i,j]
        }
        
        if(lower_marks[i,j]){
          lower_matrix[i,j]=vol_data_today[i,j]
        }
        
      }
    }
    
    upper <- round(upper_matrix*100,digits = 2)
    lower <- round(lower_matrix*100,digits = 2)
    return(list(vol_data_today,upper,lower, upper_marks,lower_marks))
  })
  pars <- reactive({
    
    SABR_Parms<-matrix(nrow=nrow( sub_vol_parms()[[1]]),ncol=3)
    SABR_Volatilities<-matrix(nrow=nrow(sub_vol_parms()[[1]]),ncol=ncol(sub_vol_parms()[[1]]))
    i <- as.numeric(input$select)
    SABR<-FitSABR(0.5,sub_vol_parms()[[2]][i,],sub_vol_parms()[[2]][i,],Ex,sub_vol_parms()[[1]][i,])
    SABR_Parms[i,]<-SABR$par
    
    return(SABR$par)
    
  })
  parsx <- reactive({
    
    SABR_Parms<-matrix(nrow=nrow( sub_vol_parmsx()[[1]]),ncol=3)
    SABR_Volatilities<-matrix(nrow=nrow(sub_vol_parmsx()[[1]]),ncol=ncol(sub_vol_parmsx()[[1]]))
    i <- as.numeric(input$selectx)
    SABR<-FitSABR(0.5,sub_vol_parmsx()[[2]][i,],sub_vol_parmsx()[[2]][i,],Ex,sub_vol_parmsx()[[1]][i,])
    SABR_Parms[i,]<-SABR$par
    
    return(SABR$par)
    
  })
  data_tomorrow <- reactive({
    location <- which(as.Date(Vols_data$Date)==as.Date(input$date_new))+1
    vols_sum<-Vols_data[location,-1]
    strikes_sum<- Strikes_data[location,-1]
    
    return(list(vols_sum,strikes_sum))
  })
  price_swaption <- reactive({
     date <- input$date_new
     data<-load_data(date,Vols_data,Strikes_data)
     vol_data_today<-t(as.data.frame(data[1]))/100
     K_data_today<-t(as.data.frame(data[2]))/100
     tomorrow_data<-data_tomorrow()
     vols_tomorrow<-tomorrow_data[[1]]/100
     strikes_tomorrow<-tomorrow_data[[2]]/100
     vols_tomorrow<-matrix(unlist(vols_tomorrow),nrow = 15, ncol = 18, byrow = FALSE )#
     strikes_tomorrow<-matrix(unlist(strikes_tomorrow),nrow = 15, ncol = 18, byrow = FALSE )#
     # if(date_tomorrow==())
     
     tom_upper_vols<-matrix(nrow=15,ncol=18)
     tom_upper_k<-matrix(nrow=15,ncol=18)
     tom_lower_vols<-matrix(nrow=15,ncol=18)
     tom_lower_k<-matrix(nrow=15,ncol=18)
     today_upper_vols<-matrix(nrow=15,ncol=18)
     today_upper_k<-matrix(nrow=15,ncol=18)
     today_lower_vols<-matrix(nrow=15,ncol=18)
     today_lower_k<-matrix(nrow=15,ncol=18)
     upper_marks<-moni_new()[[4]]
     lower_marks<-moni_new()[[5]]
     for(i in 1:nrow(vols_tomorrow)){
       for(j in 1:ncol(vols_tomorrow)){
         if(upper_marks[i,j]){
           tom_upper_vols[i,j]=vols_tomorrow[i,j]
           tom_upper_k[i,j]=strikes_tomorrow[i,j]
           today_upper_vols[i,j]=vol_data_today[i,j]
           today_upper_k[i,j]=K_data_today[i,j]
         }
         
         if(lower_marks[i,j]){
           tom_lower_vols[i,j]=vols_tomorrow[i,j]
           tom_lower_k[i,j]=strikes_tomorrow[i,j]
           today_lower_vols[i,j]=vol_data_today[i,j]
           today_lower_k[i,j]=K_data_today[i,j]
         }
         
       }
     }
     today_upper_price<-calculate_bs(today_upper_vols,today_upper_k,date,yield_data)
     today_upper_price_c<-today_upper_price
     today_lower_price<-calculate_bs(today_lower_vols,today_lower_k,date,yield_data)
     today_lower_price_c<-today_lower_price
     tomorrow_upper_price<-calculate_bs(tom_upper_vols,tom_upper_k,as.Date(date)+1,yield_data)
     tomorrow_lower_price<-calculate_bs(tom_lower_vols,tom_lower_k,as.Date(date)+1,yield_data)
    
     upper_price_sum<-today_upper_price-tomorrow_upper_price
     lower_price_sum<-tomorrow_lower_price-today_lower_price
     upper_n<-length(upper_marks[upper_marks==T])
     lower_n<-length(lower_marks[lower_marks==T])
     today_cap_sum <- today_upper_price + today_lower_price
     today_cap_sum[is.na(today_cap_sum)]=0
     upper_price_sum[is.na(upper_price_sum)]=0
     lower_price_sum[is.na(lower_price_sum)]=0
     # rate<-(sum(upper_price_sum)+sum(lower_price_sum))/sum(today_cap_sum)
     
     upper_return <- (today_upper_price-tomorrow_upper_price)/today_upper_price
     lower_return <- (tomorrow_lower_price-today_lower_price)/today_lower_price
     upper_return[is.na(upper_return)]=0
     lower_return[is.na(lower_return)]=0
     return_div <- (upper_return + lower_return)*100
     notional<-as.numeric(input$notional)
     return_div[return_div==0]=NA
     return_div <- round(return_div,digits = 2)
     today_upper_price_c[is.na(today_upper_price_c)]=0
     today_lower_price_c[is.na(today_lower_price_c)]=0
     
     
     upper_total<-sum(upper_price_sum)/sum(today_upper_price_c)
     lower_total<-sum(lower_price_sum)/sum(today_lower_price_c)
     
     upper_short<-sum(today_upper_price_c*notional)
     lower_long<-sum(today_lower_price_c*notional)
     total_notional<-as.numeric(upper_n+lower_n)*notional
     
     upper_profit<-upper_total*upper_n*notional
     lower_profit<-lower_total*lower_n*notional
     profit<-upper_profit+lower_profit
     rate<-profit/total_notional
     return(list(today_upper_price,today_lower_price,profit,rate*100,return_div,
                 upper_profit,lower_profit,upper_total*100,lower_total*100,total_notional,upper_short,lower_long))
     
   })
##### output-----  

  output$plot1 <- renderPlotly({
    
    if(weekdays(as.Date(input$date1))=="Sunday" | weekdays(as.Date(input$date1))=="Saturday"){
      
      p <-plotly_empty()%>%
        layout(title = paste0("Swaption Volatilities on last Friday,Date:",as.Date(paste(year(input$date1), 
                                                                                         week(input$date1)-1, 5, sep="-"), "%Y-%U-%u")),
               scene = list(
                 yaxis = list(title = "Tenor"),
                 xaxis = list(title = "Expiry"),
                 zaxis = list(title = "Vols")))
      if(input$flag1)
        p <- add_surface(p ,x =expiry , y = c(tenor,31,32,33) , z =  volatiltiy()[[2]], type = "surface", height = 900, width = 1200)
      if(input$flag2)    
        p <- add_surface(p, x = expiry, y =c(tenor,31,32,33),z=  volatiltiy()[[1]], type = "surface", height = 900, width = 1200)
      
      p}
    else{
    
      p <-plotly_empty()%>%
        layout(title = paste0("Swaption Volatilities on ,Date:",as.Date(input$date1)),
               scene = list(
                 yaxis = list(title = "Tenor"),
                 xaxis = list(title = "Expiry"),
                 zaxis = list(title = "Vols")))
      if(input$flag1)
        p <- add_surface(p ,x =expiry , y = c(tenor,31,32,33) , z =  volatiltiy()[[2]], type = "surface", height = 900, width = 1200)
      if(input$flag2)    
        p <- add_surface(p, x =expiry , y = c(tenor,31,32,33) , z=  volatiltiy()[[1]], type = "surface", height = 900, width = 1200)
      
      p
    
 

    }
    })
  output$plot2 <- renderPlot({
    if(as.Date(input$date2) %in% as.Date(Vols_data[,1])){
    par(mfrow=c(3,5))
    SABR_Parms<-matrix(nrow=nrow( sub_vol()[[1]]),ncol=3)
    SABR_Volatilities<-matrix(nrow=nrow(sub_vol()[[1]]),ncol=ncol(sub_vol()[[1]]))
    for(i in 1:15){
      SABR<-FitSABR(0.5,sub_vol()[[2]][i,],sub_vol()[[2]][i,],Ex,sub_vol()[[1]][i,])
      SABR_Parms[i,]<-SABR$par
      SABR_Volatilities[i,]<-SABR$estvols
      plot(x=Ex,y=sub_vol()[[1]][i,],type="l",xlab="expiry",
           ylab="vols",main=paste0("X year into ",tenor[i]," year swaption"),
           width=4, height=3)
      lines(x=Ex,y=SABR$estvols,type="l",col=2)
    }
    }
    else{
    
      par(mfrow=c(3,5))
      SABR_Parms<-matrix(nrow=nrow( sub_vol()[[1]]),ncol=3)
      SABR_Volatilities<-matrix(nrow=nrow(sub_vol()[[1]]),ncol=ncol(sub_vol()[[1]]))
      for(i in 1:15){
        SABR<-FitSABR(0.5,sub_vol()[[2]][i,],sub_vol()[[2]][i,],Ex,sub_vol()[[1]][i,])
        SABR_Parms[i,]<-SABR$par
        SABR_Volatilities[i,]<-SABR$estvols
        plot(x=Ex,y=sub_vol()[[1]][i,],type="l",xlab="expiry",
             ylab="vols",main=paste0("X year into ",tenor[i]," year swaption"),
             width=4, height=3)
        lines(x=Ex,y=SABR$estvols,type="l",col=2)
      }
      
    }

  })
  output$plot3 <- renderPlotly({
    if(as.Date(input$date3) %in% as.Date(Vols_data[,1])){
    p<-plot_ly(x =c(tenor,31,32,33) , y = expiry  , z =moni_vol()[[1]]*100, type = "surface",showscale = FALSE)%>%
      layout(showlegend = TRUE)%>%
      # add_surface(p,x=tenor,y=rep(expiry[16],15),z=as.matrix(ob_vols[,16]*100),type="surface")%>%
      add_markers(x = tenor_plus, y =rep(expiry,15)  ,z = as.numeric(split(t(moni_vol()[[3]]*100), rep(1))[[1]]),
                  mode="markers",name="Lower")%>%
      add_markers(x = tenor_plus, y =rep(expiry,15) ,z = as.numeric(split(t(moni_vol()[[2]]*100), rep(1))[[1]]),
                  mode="markers",name="Upper")%>%
      
      layout(title = paste0("Swaption Volatilities ",input$date3),
             scene = list(
               yaxis = list(title = "Expiry"),
               xaxis = list(title = "Tenor"),
               zaxis = list(title = "Vols")))
    p
    }
    else{
      p<-plotly_empty()%>%
        layout(title = paste0(input$date3," Data Not Available  "))
      p
    }
  })
  output$Alpha <- renderValueBox({
    valueBox(
      paste0(round(pars()[1],4)) , h2("Alpha"), icon = icon("cog"),
      color = "purple"
    )
  })
  output$Rho <- renderValueBox({
    valueBox(
      paste0(round(pars()[2],4)) , h2("Rho"), icon = icon("cog"),
      color = "purple"
    )
  })
  output$Nu <- renderValueBox({
    valueBox(
      paste0(round(pars()[3],4) ), h2("Nu"), icon = icon("cog"),
      color = "purple"
    )
  })
  output$profit <- renderValueBox({
    valueBox(
      paste0(round(price_swaption()[[3]],digits=2)) , "Profit", icon = icon("cog"),
      color = "purple",width=4
    )
  })
  # list(today_upper_price,today_lower_price,profit,rate,return_div,upper_profit,lower_profit,upper_total,lower_total)
  output$rate <- renderValueBox({
    valueBox(
      paste0(round(price_swaption()[[4]],digits = 2),"%") , h2("Total Return"), icon = icon("percent"),
      color = "purple",width=1
    )
  })
  output$red <- renderValueBox({
    valueBox(
      paste0(round(price_swaption()[[8]],digits = 2),"%") , h2("Red Return"), icon = icon("percent"),
      color = "red",width=1
    )
  })
  output$profit <- renderValueBox({
    valueBox(
      paste0(round(price_swaption()[[3]],digits = 2)) , h2("Total Porfit"), icon = icon("money"),
      color = "purple",width=1
    )
  })
  output$r_profit <- renderValueBox({
    valueBox(
      paste0(round(price_swaption()[[6]],digits = 2)) , h2("Red Profit"), icon = icon("money"),
      color = "red",width=1
    )
  })
  output$g_profit <- renderValueBox({
    valueBox(
      paste0(round(price_swaption()[[7]],digits = 2)) , h2("Green Profit"), icon = icon("money"),
      color = "green",width=1
    )
  })
  output$green <- renderValueBox({
    valueBox(
      paste0(round(price_swaption()[[9]],digits = 2),"%") , h2("Green Return"), icon = icon("percent"),
      color = "green",width=1
    )
  })
  output$long_amount <- renderValueBox({
    valueBox(
      paste0(round(price_swaption()[[12]],digits=2)) , h2("Long"), icon = icon("money"),
      color = "yellow"
    )
  })
  output$short_amount <- renderValueBox({
    valueBox(
      paste0(round(price_swaption()[[11]],digits=2)) , h2("Short"), icon = icon("money"),
      color = "yellow"
    )
  })
  output$total_amount <- renderValueBox({
    valueBox(
      paste0(round(price_swaption()[[12]]-price_swaption()[[11]],digits=2)) , h2("Total"), icon = icon("money"),
      color = "yellow"
    )
  })
  output$total_notional <- renderValueBox({
    valueBox(
      paste0(round(price_swaption()[[10]],digits=2)) , h2("Total notional"), icon = icon("bank"),
      color = "yellow"
    )
  })
  # list(today_upper_price,today_lower_price,profit,rate*100,return_div,
       # upper_profit,lower_profit,upper_total*100,lower_total*100,total_notional,upper_short,lower_long)
  output$plotgg <- renderPlot({
    ggplot(  moni_vol()[[4]], aes(x=df_Ex, y=df_tenor)) + 
      geom_tile( fill= "white", color = "white") +
      geom_text(aes(label=df_upperlower, color=  df_color<1.5)) +
      scale_color_manual( name= "Legend",labels=c("Upper","Lower"),values=c("red", "green")) +
      # scale_fill_manual(labels = c("FALSE" = "Sell", "TRUE" = "Buy"), 
      #                   values = c('red', 'green')) + 
      
      scale_y_discrete(name ="Tenor", 
                       limits=c("1","2","3","4","5","6","7","8","9","10","12","15","20","25","30"),
                       labels=c("1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y","12Y","15Y","20Y","25Y","30Y"))+
      scale_x_discrete(name="Expriy",breaks=as.factor(c(1/12,3/12,6/12,9/12,1,2,3,4,5,6,7,8,9,10,15,20,25,30)),
                       labels=c("1M","3M","6M","9M","1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y","15Y","20Y","25Y","30Y"))+
      theme(axis.text = element_text(face = "bold")) +
      ggtitle(
        paste0(
        "Trading suggestions on ",input$date3)
        )
  })
  output$plot_new <- renderPlotly({
    #list(vol_data_today,upper,lower, upper_marks,lower_marks)
   
    
    p<-plot_ly(x =  expiry, y =c(tenor,31,32,33)  , z =moni_new()[[1]]*100, type = "surface",showscale = FALSE)%>%
      layout(showlegend = TRUE)%>%
      # add_surface(p,x=tenor,y=rep(expiry[16],15),z=as.matrix(ob_vols[,16]*100),type="surface")%>%
      add_markers(x = rep(expiry,15), y = tenor_plus ,z = as.numeric(split(t(moni_new()[[2]]), rep(1))[[1]]),
                  mode="markers",name="Upper")%>%
      add_markers(x = rep(expiry,15), y = tenor_plus ,z = as.numeric(split(t(moni_new()[[3]]), rep(1))[[1]]),
                  mode="markers",name="Lower")%>%
      
      layout(title = paste0("Swaption Volatilities ",input$date_new),
             scene = list(
               xaxis = list(title = "Expiry"),
               yaxis = list(title = "Tenor"),
               zaxis = list(title = "Vols")))
    p
  })
  output$plot4 <- renderPlotly({
    
    if(as.Date(input$date2) %in% as.Date(Vols_data[,1])){
      
      SABR_Parms<-matrix(nrow=nrow( sub_vol_parms()[[1]]),ncol=3)
      SABR_Volatilities<-matrix(nrow=nrow(sub_vol_parms()[[1]]),ncol=ncol(sub_vol_parms()[[1]]))
      i <- as.numeric(input$select)
        SABR<-FitSABR(0.5,sub_vol_parms()[[2]][i,],sub_vol_parms()[[2]][i,],Ex,sub_vol_parms()[[1]][i,])
        SABR_Parms[i,]<-SABR$par
        SABR_Volatilities[i,]<-SABR$estvols
      #   # plot(x=Ex,y=sub_vol()[[1]][i,],type="l",xlab="expiry",
      #   #      ylab="vols",main=paste0("X year into ",tenor[i]," year swaption"),
      #   #      width=4, height=3)
      #   # lines(x=Ex,y=SABR$estvols,type="l",col=2)
        p <- plot_ly( x = Ex, y = sub_vol_parms()[[1]][i,], name = 'Observed', type = 'scatter', mode = 'lines') %>%
          add_trace(y = SABR$estvols, name = 'SABR model', mode = 'lines')%>%
          layout(
            title = paste0("X in to ",tenor[as.numeric(input$select)],"Y Swaption Volatilities on ,Date:",as.Date(input$date2)),
                 scene = list(
                   xaxis = list(title = "Expiry"),
                   yaxis = list(title = "Vols")))
      
    }
    else{
      
      
      SABR_Parms<-matrix(nrow=nrow( sub_vol_parms()[[1]]),ncol=3)
      SABR_Volatilities<-matrix(nrow=nrow(sub_vol_parms()[[1]]),ncol=ncol(sub_vol_parms()[[1]]))
      i <- as.numeric(input$select)
        SABR<-FitSABR(0.5,sub_vol_parms()[[2]][i,],sub_vol_parms()[[2]][i,],Ex,sub_vol_parms()[[1]][i,])
        SABR_Parms[i,]<-SABR$par
        SABR_Volatilities[i,]<-SABR$estvols
        p <- plot_ly( x = Ex, y = sub_vol_parms()[[1]][i,], name = 'Observed', type = 'scatter', mode = 'lines') %>%
          add_trace(y = SABR$estvols, name = 'SABR model', mode = 'lines')%>%
        layout(
          title = paste0("X in to ",tenor[as.numeric(input$select)],"Y Swaption Volatilities on last Friday,Date:",as.Date(paste(year(input$date2), 
                                                                                    week(input$date2)-1, 5, sep="-"), "%Y-%U-%u")),
            scene = list(
              xaxis = list(title = "Expiry"),
              yaxis = list(title = "Vols")))
      #   
      # 
        
               
    }
    
  })
  output$plotgg_new <- renderPlot({
    EX_char <- c(1/12,3/12,6/12,9/12,1,2,3,4,5,6,7,8,9,10,15,20,25,30)
    EX_char <- as.factor(EX_char)
    Ex_char_plus <- rep(EX_char,15)
    tenor_char <-c(1,2,3,4,5,6,7,8,9,10,12,15,20,25,30)
    tenor_char <- as.factor(tenor_char)
    tenor_char_plus<-c()
    for(i in tenor_char){
      tenor_char_plus<-c(tenor_char_plus,rep(i,18))
    }
    upper_plot<-split(t(moni_new()[[2]]), rep(1))[[1]]
    lower_plot<-split(t(moni_new()[[3]]), rep(1))[[1]]
    
    upper_lower<-c(upper_plot,lower_plot)
    upper_lower<-as.data.frame(upper_lower)
    Ex_char_plus<-as.data.frame(rep(Ex_char_plus,2))
    
    tenor_char_plus<-as.data.frame(rep(tenor_char_plus,2))
    col_flag <- as.data.frame(c(rep(1,270),rep(2,270)))
    upper_lower_plot <- cbind(upper_lower,Ex_char_plus,tenor_char_plus,col_flag)
    upper_lower_plot <- as.data.frame(upper_lower_plot)
    colnames(upper_lower_plot) <- c("df_upperlower","df_Ex","df_tenor","df_color")
    
    #ggplot---------
    
    ggplot(  upper_lower_plot, aes(x=df_Ex, y=df_tenor)) + 
      geom_tile( fill= "white", color = "white") +
      geom_text(aes(label=df_upperlower, color=  df_color<1.5)) +
      scale_color_manual( name= "Legend",labels=c("Lower","Upper"),values=c("green", "red")) +
      # scale_fill_manual(labels = c("FALSE" = "Sell", "TRUE" = "Buy"), 
      #                   values = c('red', 'green')) + 
      
      scale_y_discrete(name ="Tenor", 
                       limits=c("1","2","3","4","5","6","7","8","9","10","12","15","20","25","30"),
                       labels=c("1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y","12Y","15Y","20Y","25Y","30Y"))+
      scale_x_discrete(name="Expriy",breaks=as.factor(c(1/12,3/12,6/12,9/12,1,2,3,4,5,6,7,8,9,10,15,20,25,30)),
                       labels=c("1M","3M","6M","9M","1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y","15Y","20Y","25Y","30Y"))+
      theme(axis.text = element_text(face = "bold")) +
      ggtitle(
        paste0(
          "Trading suggestions on ",input$date_new," With ",input$c*100,"% Confidence interval ")
      )
    
    
  })
  output$plotgg_price <- renderPlot({
    
    EX_char <- c(1/12,3/12,6/12,9/12,1,2,3,4,5,6,7,8,9,10,15,20,25,30)
    EX_char <- as.factor(EX_char)
    Ex_char_plus <- rep(EX_char,15)
    tenor_char <-c(1,2,3,4,5,6,7,8,9,10,12,15,20,25,30)
    tenor_char <- as.factor(tenor_char)
    tenor_char_plus<-c()
    for(i in tenor_char){
      tenor_char_plus<-c(tenor_char_plus,rep(i,18))
    }
    today_upper_price<-price_swaption()[[1]]
    today_lower_price<-price_swaption()[[2]]
    upper_plot<-split(t(round(today_upper_price,digits=4)), rep(1))[[1]]
    lower_plot<-split(t(round(today_lower_price,digits=4)), rep(1))[[1]]
    upper_lower<-c(upper_plot,lower_plot)
    upper_lower<-as.data.frame(upper_lower)
    Ex_char_plus<-as.data.frame(rep(Ex_char_plus,2))
    
    tenor_char_plus<-as.data.frame(rep(tenor_char_plus,2))
    col_flag <- as.data.frame(c(rep(1,270),rep(2,270)))
    upper_lower_plot <- cbind(upper_lower,Ex_char_plus,tenor_char_plus,col_flag)
    upper_lower_plot <- as.data.frame(upper_lower_plot)
    colnames(upper_lower_plot) <- c("df_upperlower","df_Ex","df_tenor","df_color")
    
    #ggplot---------
    
    ggplot(  upper_lower_plot, aes(x=df_Ex, y=df_tenor)) + 
      geom_tile( fill= "white", color = "white") +
      geom_text(aes(label=df_upperlower, color=  df_color<1.5),size=4) +
      scale_color_manual( name= "Legend",labels=c("Lower","Upper"),values=c("green", "red")) +
      # scale_fill_manual(labels = c("FALSE" = "Sell", "TRUE" = "Buy"), 
      #                   values = c('red', 'green')) + 
      
      scale_y_discrete(name ="Tenor", 
                       limits=c("1","2","3","4","5","6","7","8","9","10","12","15","20","25","30"),
                       labels=c("1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y","12Y","15Y","20Y","25Y","30Y"))+
      scale_x_discrete(name="Expriy",breaks=as.factor(c(1/12,3/12,6/12,9/12,1,2,3,4,5,6,7,8,9,10,15,20,25,30)),
                       labels=c("1M","3M","6M","9M","1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y","15Y","20Y","25Y","30Y"))+
      theme(axis.text = element_text(face = "bold")) +
      ggtitle(
        paste0(
          "Upper and Lower Swaption price on Date: ",input$date_new," With ",input$c*100,"% Confidence interval ")
      )+
      theme(text = element_text(size=15))+
      theme(plot.title = element_text(size=18))+
      theme(legend.text=element_text(size=12))+
      theme(axis.title = element_text(size=15))
    
    


    
    
    
    
    
  })
  output$plotx <- renderPlotly({
    
    if(as.Date(input$datex) %in% as.Date(Vols_data[,1])){
      
      SABR_Parms<-matrix(nrow=nrow( sub_vol_parmsx()[[1]]),ncol=3)
      SABR_Volatilities<-matrix(nrow=nrow(sub_vol_parmsx()[[1]]),ncol=ncol(sub_vol_parmsx()[[1]]))
      i <- as.numeric(input$selectx)
      # SABR<-FitSABR(input$betax,sub_vol_parmsx()[[1]],input$rhox,Ex,input$nux)
      # SABR_Parms[i,]<-SABR$par
      SABR_Volatilities[i,]<-SABRVol(input$alphax,input$betax,input$rhox,input$nux,Ex,sub_vol_parmsx()[[2]][i,],sub_vol_parmsx()[[2]][i,])
      #   # plot(x=Ex,y=sub_vol()[[1]][i,],type="l",xlab="expiry",
      #   #      ylab="vols",main=paste0("X year into ",tenor[i]," year swaption"),
      #   #      width=4, height=3)
      #   # lines(x=Ex,y=SABR$estvols,type="l",col=2)
      p <- plot_ly( x = Ex, y = sub_vol_parmsx()[[1]][i,], name = 'Observed', type = 'scatter', mode = 'lines') %>%
        add_trace(y = SABR_Volatilities[i,], name = 'SABR model', mode = 'lines')%>%
        layout(
          title = paste0("X in to ",tenor[as.numeric(input$selectx)],"Y Swaption Volatilities on ,Date:",as.Date(input$datex)),
          scene = list(
            xaxis = list(title = "Expiry"),
            yaxis = list(title = "Vols")))
      
    }
    else{
      
      
      SABR_Parms<-matrix(nrow=nrow( sub_vol_parmsx()[[1]]),ncol=3)
      SABR_Volatilities<-matrix(nrow=nrow(sub_vol_parmsx()[[1]]),ncol=ncol(sub_vol_parmsx()[[1]]))
      i <- as.numeric(input$selectx)
      # SABR<-FitSABR(0.5,sub_vol_parmsx()[[2]][i,],sub_vol_parmsx()[[2]][i,],Ex,sub_vol_parmsx()[[1]][i,])
      # SABR_Parms[i,]<-SABR$par
      SABR_Volatilities[i,]<-SABRvol(input$alphax,input$betax,input$rhox,input$nux,Ex,sub_vol_parmsx()[[2]][i,],sub_vol_parmsx()[[2]][i,])
      p <- plot_ly( x = Ex, y = sub_vol_parmsx()[[1]][i,], name = 'Observed', type = 'scatter', mode = 'lines') %>%
        add_trace(y = SABR_Volatilities[i,], name = 'SABR model', mode = 'lines')%>%
        layout(
          title = paste0("X in to ",tenor[as.numeric(input$selectx)],"Y Swaption Volatilities on last Friday,Date:",as.Date(paste(year(input$datex), 
                                                                                                                                 week(input$datex)-1, 5, sep="-"), "%Y-%U-%u")),
          scene = list(
            xaxis = list(title = "Expiry"),
            yaxis = list(title = "Vols")))
      #   
      # 
      
      
    }
    
  })
  output$plotgg_return <- renderPlot({
   
    test_matrix <- price_swaption()[[5]]
    # test_matrix <- matrix(c(1:270),ncol=18,nrow=15)
    low_limit <- -input$return1
    high_limit <- input$return2
    test_list<-split(t(test_matrix), rep(1))[[1]]
    test_flag <- unlist(lapply(test_list,function(x){
      if(!is.na(x)&&x<low_limit) return(1)
      else if (!is.na(x)&&low_limit<=x && x<high_limit) return(2)
      else if(!is.na(x)&&x>=high_limit)return(3)
      else return(NA)
      
    }))
    
    
    EX_char <- round(c(1/12,3/12,6/12,9/12,1,2,3,4,5,6,7,8,9,10,15,20,25,30),digits=3)
    EX_char <- as.factor(EX_char)
    Ex_char_plus <- rep(EX_char,15)
    tenor_char <-c(1,2,3,4,5,6,7,8,9,10,12,15,20,25,30)
    tenor_char <- as.factor(tenor_char)
    tenor_char_plus<-c()
    for(i in tenor_char){
      tenor_char_plus<-c(tenor_char_plus,rep(i,18))
    }
    
    df_gg <- data.frame(test_list,Ex_char_plus,tenor_char_plus,test_flag)
    colnames(df_gg) <- c("df_upperlower","df_Ex","df_tenor","df_color")
    
    ggplot(  df_gg, aes(x=df_Ex, y=df_tenor)) +
      geom_tile( fill= "white", color = "white") +
      geom_text(aes(label=df_upperlower, color=as.factor(df_color)
      ),size=4.5) +
      scale_color_manual( name= "Legend",labels=c("Loss","mid","Gain"),values=c("#ff3300","black", "#0000ff")) +
      # scale_fill_manual(labels = c("Lower" = "Sell", "Upper" = "Buy"),
      #                   values = c('red', 'green')) +
      
      scale_y_discrete(name ="Tenor",
                       limits=c("1","2","3","4","5","6","7","8","9","10","12","15","20","25","30"),
                       labels=c("1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y","12Y","15Y","20Y","25Y","30Y"))+
      scale_x_discrete(name="Expriy",breaks=as.factor(c(1/12,3/12,6/12,9/12,1,2,3,4,5,6,7,8,9,10,15,20,25,30)),
                       labels=c("1M","3M","6M","9M","1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y","15Y","20Y","25Y","30Y"))+
      theme(axis.text = element_text(size=15)) +
      ggtitle(
        paste0(
          "Upper and Lower Swaption return(%) on Date: ",input$date_new," With ",input$c*100,"% Confidence interval ")
       
      )+
      theme(text = element_text(size=15))+
     theme(plot.title = element_text(size=18))+
     theme(legend.text=element_text(size=12))+
      theme(axis.title = element_text(size=15))
    
  })
  output$yield_c <- renderPlotly({
    t<-seq(from = 6 , to = 360 ,by = 6)
    Time<-c()
    for(i in 1:length(t)){
      Time[i]<-paste0(t[i],"M")
    }
    # Time <- as.factor(unlist(lapply(x,function(x) return(paste0(x,"M")))))
    date<-as.Date(input$date_new)
    if(weekdays(date)=="Friday"){date_to<-date+3}
    else{date_to<-date+1}
    # date<-as.Date("2016-7-14")
    discount_factor<-as.numeric(yield_curve(date,yield_data))
    discount_factor_to<-as.numeric(yield_curve(date+1,yield_data))
    p<-plot_ly(x=t,y=discount_factor,type = 'scatter',name=date, mode = 'lines')%>%
      add_trace(x=t,y=discount_factor_to,name=date_to)%>%
      layout(
        title = paste0("Yield Curve on Date: ",date),
        
        scene = list(
          xaxis = list(title = "Month"),
          yaxis = list(title = "Discount factor")))
    p
    
})
}
shinyApp(ui, server)