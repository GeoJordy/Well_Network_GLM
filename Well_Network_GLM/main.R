rm(list=objects())
memory.limit(size = 16273*10)
# Load up data
{
  this_fp <- function(){return(rstudioapi::getActiveDocumentContext()$path)}
  this_dir <- function(){return(dirname(this_fp()))}
  root <- this_dir()
  setwd(root)
  load("data_in/mk_data.Rdata")
}
# load up dependencies
{
  library(jagsUI)
  library(jagshelper)
  library(dplyr)
  get_month <- function(date){
    as.numeric(as.character(date, format="%m"))
  }
}
# stratify by depth and prepare data for use in Jags model
{
  tab_tce <- tab[tab$D1 > 0,]
  plot(tab_tce$D1, tab_tce$NormalizedResult)
  prepare_jags_data <- function(tab, analyte){
    tab <- tab[tab$Analyte == analyte,]
    data_rel <- data.frame(
      Date = tab$`Date Sampled`,
      Well = tab$Location,
      Depth = tab$D1,
      Z = tab$NormalizedResult,
      Det = tab$is_detected,
      DL = tab$DetectionLimit,
      iwells = NA,
      idepth = NA
    )
    data_rel <- dplyr::arrange(data_rel, Date)
    start_date <- data_rel$Date[1]
    start_month <- get_month(start_date)
    data_rel$Day <- as.numeric(
      floor(difftime(data_rel$Date, start_date, units = "days"))
    ) + 1
    data_rel$Month <- get_month(data_rel$Date)
    wells <- sort(unique(data_rel$Well))
    for (i in 1:length(wells)){
      well <- wells[i]
      data_rel$iwells[data_rel$Well==well] <- i 
    }
    print(paste(sum(data_rel$Depth < 0), "many measurements with depth < 0"))
    print("removing them from the data set")
    data_rel <- data_rel[data_rel$Depth > 0,]
    print(paste(sum(data_rel$DL < 0), "many measurements without detection limits"))
    print("removing them from the data set")
    data_rel <- data_rel[data_rel$DL > 0,]
    data_rel$idepth[data_rel$Depth < 25] <- 1
    data_rel$idepth[data_rel$Depth >= 25 & data_rel$Depth < 50] <- 2
    data_rel$idepth[data_rel$Depth >= 50 & data_rel$Depth < 75] <- 3
    data_rel$idepth[data_rel$Depth >= 75 & data_rel$Depth < 100] <- 4
    data_rel$idepth[data_rel$Depth >= 100 & data_rel$Depth < 125] <- 5
    data_rel$idepth[data_rel$Depth >= 125] <- 6
    jags_data <- list(
      n_depths = length(unique(data_rel$idepth)), 
      n_wells = length(unique(data_rel$iwells)),
      n_d = sum(data_rel$Det==1),
      Y_d = data_rel$Z[data_rel$Det==1],
      Time_d = data_rel$Day[data_rel$Det==1],
      Depth_d = data_rel$idepth[data_rel$Det==1],
      Well_d = data_rel$iwells[data_rel$Det==1],
      n_u = sum(data_rel$Det==0),
      DL_u = data.frame(0.01, data_rel$DL[data_rel$Det==0]),
      Time_u = data_rel$Day[data_rel$Det==0],
      Depth_u = data_rel$idepth[data_rel$Det==0],
      Well_u = data_rel$iwells[data_rel$Det==0],
      ones = rep(1, sum(data_rel$Det==0))
    )
    print(paste("Time NAs:",  sum(is.na(jags_data$Time))))
    print(paste("Month NAs:" , sum(is.na(jags_data$Month))))
    print(paste("Depth NAs" , sum(is.na(jags_data$Depth))))
    print(paste("Well NAs" , sum(is.na(jags_data$Well))))
    print(paste("DL NAs" , sum(is.na(jags_data$DL))))
    print(paste("Z NAs" , sum(is.na(jags_data$Z))))
    return(jags_data)
  }
  jags_data_tce <- prepare_jags_data(tab, "Trichloroethene (TCE)")
}
# build jags model
{
  cat(
    'model{
      # Likelihood -- for observed data
      for (i in 1:n_d){
        mu_d[i] <- beta1[Depth_d[i]] + beta2[Well_d[i]] + beta3[Well_d[i]]*Time_d[i]
        Y_d[i] ~ dlnorm(
          mu_d[i],
          tauw[Well_d[i]]
        )
        # res_d[i] <- log(Y_d[i]) - mu_d[i] 
      }
      # Likelihood -- for unobserved data (i.e. non-detects)
      for (i in 1:n_u){
        ones[i] ~ dinterval(Z[i], DL_u[i,])
        mu_u[i] <- beta1[Depth_u[i]] + beta2[Well_u[i]] + beta3[Well_u[i]]*Time_u[i]
        Z[i] ~ dlnorm(
          mu_u[i],
          tauw[Well_u[i]]
        )
        # res_u[i] <- log(Z[i]) - mu_u[i]
      }
      for (i in 1:n_depths){
        beta1[i] ~ dnorm(0, 0.001)
      }
      for (i in 1:n_wells){
        beta2[i] ~ dnorm(0, 0.001)
        beta3[i] ~ dnorm(0, 0.001)
        tauw[i] <- pow(1/sigw[i], 2)
        sigw[i] ~ dexp(0.001)
      }
    }', 
    file='jags_model'
  )
}
# run jags models
{
  params_out <- c(
    "beta0",
    "beta1",
    "beta2",
    "beta3",
    "sigw",
    "mu_d",
    "mu_u",
    "Z"
    # "res_u",
    # "res_d"
  )
  niter <- 500000
  ncores <- 6
}
{
  start_time <- Sys.time()
  jags_out_tce <- jagsUI::jags(model.file='jags_model', 
                               data=jags_data_tce,
                               parameters.to.save=params_out,
                               n.chains=ncores, 
                               parallel=T, 
                               n.iter=niter,
                               n.burnin=niter/2, 
                               n.thin=niter/2000)
  save(jags_out_tce, file="jags_out/jags_out_tce.Rdata")
  print(Sys.time() - start_time)
}  