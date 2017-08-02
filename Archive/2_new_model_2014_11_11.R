##################################################################################
source("F:/_SM/SMAO/Felles/autostat/Libraries/r_tools/r_tools.R")
req(c("sfsmisc","deSolve","data.table","boot","rmarkdown","knitr"))

setwd("F:/_SM/_SM-Felles/Utbrudd/Ebola 2014/Ebola_models_updated")

days.break <- 28*2
time.incubation=11.4
time.to.report=5
time.infectedL=5+11.8
time.infectedD=5+4.2
cfr=0.7
reported=0.5


today.character <- as.character(as.Date(format(Sys.time(), "%Y-%m-%d")))
try(dir.create(today.character),TRUE)

##################################################################################
##################################################################################
# Let's set up some initial conditions at time t=0
##################################################################################
proper <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}

calculate.rates <- function(t, x, vparameters, reported,beta.replacement=-99,time.to.report=5){
  S = x[1]  # the value of S at time t
  E = x[2]
  I = x[3]  # the value of I at time t
  R = x[4]  # the value of R at time t
  D = x[5]
  T = x[6]
  TR = x[7]
  TD = x[8]
  
  R0=vparameters["R0"]
    
  mu = 0
  a = 1/time.incubation         # approximate incubation period of influenza in days^{-1} 
  
  time.I.R = time.infectedL
  time.I.D = time.infectedD
  time.T.TR = time.infectedL-time.to.report
  time.T.TD = time.infectedD-time.to.report
  time.I.T = time.to.report
  
  gamma.I.R = 1/time.infectedL        # approximate average recovery period when living in days^{-1} 
  gamma.I.D = 1/time.infectedD        # approximate average recovery period when dead in days^{-1} 
  gamma.T.TR = 1/(time.infectedL-time.to.report)        # approximate average recovery period when living in days^{-1} 
  gamma.T.TD = 1/(time.infectedD-time.to.report)        # approximate average recovery period when dead in days^{-1} 
  gamma.I.T = 1/(time.to.report)
  
  gamma.I.R = gamma.I.R*((1-reported)*(1-cfr))
  gamma.I.D = gamma.I.D*((1-reported)*cfr)
  gamma.T.TR = gamma.T.TR*(reported*(1-cfr))
  gamma.T.TD = gamma.T.TD*(reported*cfr)
  gamma.I.T = gamma.I.T*reported
  
  time.I = time.I.T*reported + (1-reported)*((1-cfr)*time.I.R + cfr*time.I.D)
  gamma.leaving.I = 1/time.I
  gamma.I.T = gamma.leaving.I*(time.I.T*reported/time.I)
  gamma.I.R = gamma.leaving.I*((1-reported)*(1-cfr)*time.I.R/time.I)
  gamma.I.D = gamma.leaving.I*((1-reported)*(cfr)*time.I.D/time.I)
  

  time.T = ((1-cfr)*time.T.TR + cfr*time.T.TD)
  gamma.leaving.T = 1/time.T
  gamma.T.TR = gamma.leaving.T*((1-cfr)*time.T.TR/time.I)
  gamma.T.TD = gamma.leaving.T*(cfr*time.T.TD/time.I)
  
  
  if(beta.replacement==-99){
    beta  = R0*(gamma.leaving.I+mu)*(a+mu)/a   # transmission rate of influenza is estimated by solving R0=beta/gamma for beta
  } else beta=beta.replacement
  
  npop = S+E+I+R+T+TR   # the population size is always S+I+R because there are no births or deaths in the model
  dS = npop*mu - mu*S - beta*I/npop*S   # the derivative of S wrt time
  dE = beta*I/npop*S - (mu + a)*E  # the derivative of I wrt time
  
  dI = a*E - (gamma.I.R + gamma.I.D + gamma.I.T + mu)*I
  dR = gamma.I.R*I - mu*R                  # the derivative of R wrt time
  dD = gamma.I.D*I                 # the derivative of R wrt time
  
  dT = gamma.I.T*I - (gamma.T.TR+gamma.T.TD+mu)*T
  dTR =  gamma.T.TR*T - mu*TR                  # the derivative of R wrt time
  dTD = gamma.T.TD*T                  # the derivative of R wrt time
  
  out = c(gamma.leaving.I,beta,dS,dE,dI,dR,dD,dT,dTR,dTD)
  return(out)
}

SEIR.derivatives.function=function(t, x, vparameters, mod,referenceday=1,reduced.day=1,reduced.scenario=FALSE,time.to.report=5){
  for(i in 1:8) if(x[i]<0) x[i]=0

  beta.replacement = -99
  f <- data.frame(time=t)
  reported <- inv.logit(predict(mod,f))
  
  if(reduced.scenario & t>=reduced.day){
    f <- data.frame(time=referenceday)
    reported <- inv.logit(predict(mod,f))
    out = calculate.rates(t, x, vparameters, reported,time.to.report=time.to.report)
    original.gamma.leaving.I=out[1]
    original.beta=out[2]
    
    out = calculate.rates(t, x, vparameters, 0.7,time.to.report=time.to.report)
    new.gamma.leaving.I=out[1]
    
    beta.replacement = original.beta*(original.gamma.leaving.I/new.gamma.leaving.I)
    reported <- 0.7
    
  }
    
  out = calculate.rates(t, x, vparameters, reported, beta.replacement,time.to.report=time.to.report)
  out = out[-c(1:2)]
  
  list(out)
}

SEIR <- function(simulation.start,simulation.end,R0,time.incubation,time.infectedL,time.infectedD,cfr,
                 E_0=1,I_0=1,R_0=0,D_0=0,T_0=0,TR_0=0,TD_0=0,mod,reported.to.date,referenceday=1,reduced.day=1,reduced.scenario=FALSE,time.to.report=5){
  npop = 5000000
  S_0 = npop-E_0-I_0-R_0-D_0-T_0-TR_0-TD_0
  mu    = 0         # birth rate
  
  vt = seq(simulation.start,simulation.end,1)  # let's determine the values of S,I and R at times in vt
  
  vparameters = c(R0=R0)
  
  inits = c(S=S_0,E=E_0,I=I_0,R=R_0,D=D_0,T=T_0,TR=TR_0,TD=TD_0)
  sirmodel = as.data.frame(lsoda(inits, vt, SEIR.derivatives.function, vparameters,mod=mod,referenceday=referenceday,reduced.day=reduced.day,reduced.scenario=reduced.scenario,time.to.report=time.to.report))
  
  sirmodel$cases.estimated <- sirmodel$I+sirmodel$R+sirmodel$D+sirmodel$T+sirmodel$TR+sirmodel$TD
  sirmodel$cases.reported <- sirmodel$T+sirmodel$TR+sirmodel$TD
  
  return(sirmodel)
}


SEIR.against.real.data.single <- function(simulation.start,simulation.end,opts,data,mod,reported.to.date){
  predictions <- list()
  opts$dif <- 0
  opts$obs.end <- 0
  for(i in 1:nrow(opts)){
    R0=opts$R0[i]
    E_0=opts$E_0[i]
    I_0=opts$I_0[i]
    R_0=opts$R_0[i]
    D_0=opts$D_0[i]
    T_0=opts$T_0[i]
    TR_0=opts$TR_0[i]
    TD_0=opts$TD_0[i]
    
    failed <- TRUE
    try({
      sirmodel <- SEIR(simulation.start,simulation.end,
                       R0=R0,
                       time.incubation=time.incubation,
                       time.infectedL=time.infectedL,
                       time.infectedD=time.infectedD,
                       cfr=cfr,
                       E_0=E_0,
                       I_0=I_0,
                       R_0=R_0,
                       D_0=D_0,
                       T_0=T_0,
                       TR_0=TR_0,
                       TD_0=TD_0,
                       mod,reported.to.date)
      failed <- FALSE
    },FALSE)
    if(failed){
    	print("FAILED")
      opts$dif[i] <- -99
      next
    }
    
    names(sirmodel)[1] <- "day"
    res <- merge(sirmodel,data,all=TRUE)
    res$Date <- res$Date[1]+res$day-1
    res$days.to.today <- as.numeric(as.Date(format(Sys.time(), "%Y-%m-%d"))-res$Date)
    res$parameter.set <- i
    predictions[[i]] <- res[!is.na(res$cases.estimated),]
    
    testing <- res[-1,]
    testing <- testing[!is.na(testing$cases.estimated) & !is.na(testing$cases.observed),]
    
    opts$dif[i] <- mean(((testing$cases.observed-testing$cases.reported)^2))
    opts$obs.end[i] <- testing$cases.observed[nrow(testing)]
  }
  
  opts$parameter.set <- 1:nrow(opts)
  opts$d.start <- simulation.start
  opts$d.end <- simulation.end
  
  opts <- opts[opts$dif>=0,]
  
  retval <- list()
  retval[["opts"]] <- opts
  retval[["predictions"]] <- predictions
  return(retval)
}


SEIR.against.real.data.multiple <- function(simulation.start,simulation.end,opts,data,mod,reported.to.date){
  retval <- SEIR.against.real.data.single(simulation.start,simulation.end,opts,data,mod,reported.to.date)

  returned.opts <- retval[["opts"]]
  predictions <- retval[["predictions"]]
  
  returned.opts$S_end <- 0
  returned.opts$E_end <- 0
  returned.opts$I_end <- 0
  returned.opts$R_end <- 0
  returned.opts$D_end <- 0
  returned.opts$T_end <- 0
  returned.opts$TR_end <- 0
  returned.opts$TD_end <- 0
  returned.opts$time <- 0
  for(i in 1:nrow(returned.opts)){
    p <- returned.opts$parameter.set[i]
    p <- predictions[[p]]
    p <- p[nrow(p),]
    returned.opts$S_end[i] <- p$S
    returned.opts$E_end[i] <- p$E
    returned.opts$I_end[i] <- p$I
    returned.opts$R_end[i] <- p$R
    returned.opts$D_end[i] <- p$D
    returned.opts$T_end[i] <- p$T
    returned.opts$TR_end[i] <- p$TR
    returned.opts$TD_end[i] <- p$TD
    returned.opts$reported.end[i] <- p$cases.reported
    returned.opts$time[i] <- p$day
  }
  
  final.opts <- list()
  for(i in 1:nrow(returned.opts)){
    final.opts[[length(final.opts)+1]] <- returned.opts[i,]
  }
  return(final.opts)
}

calc.prob <- function(final.opts){
#  for(val in c(0.05,0.1,0.2,0.5)){
  for(val in c(0.2,0.5)){
    collapsed.final.opts <- rbindlist(final.opts)
    collapsed.final.opts <- collapsed.final.opts[,
                                                 list(dif=mean(dif)),
                                                 by=list(index)]
    ldpd <- min(collapsed.final.opts$dif)
    sigma <- sqrt(ldpd*val)
    collapsed.final.opts$dif <- collapsed.final.opts$dif/(2*(sigma^2))
    collapsed.final.opts$dif <- round(collapsed.final.opts$dif-min(collapsed.final.opts$dif),5)
    collapsed.final.opts$p <- exp(-collapsed.final.opts$dif)
    collapsed.final.opts$p <- round(collapsed.final.opts$p/sum(collapsed.final.opts$p),5)
    if(sum(collapsed.final.opts$p>0)>1) break
  }
	print("value selected:")
	print(val)
  
  #collapsed.final.opts$p[collapsed.final.opts$p!=max(collapsed.final.opts$p)] <- 0
  
  return(collapsed.final.opts)
}

reduce.options <- function(final.opts){
  for(i in 1:length(final.opts)) final.opts[[i]]$index <- i
  collapsed.final.opts <- calc.prob(final.opts)
  
  for(i in length(final.opts):1){
    final.opts[[i]] <- final.opts[[i]][,-which(names(final.opts[[i]])=="index")]
    if(collapsed.final.opts$p[i]==0){
      final.opts[[i]] <- NULL
    }
  }
  return(final.opts)
}

SEIR.second.step.against.real.data.multiple <- function(simulation.start,simulation.end,last.step,R0,data,mod){ 
  save.opts <- list()
  for(i in 1:length(last.step)){
    opts.of.interest <- last.step[[i]]
    opts.of.interest <- opts.of.interest[nrow(opts.of.interest),]
    opts <- data.frame(R0,
                       opts.of.interest$E_end,
                       opts.of.interest$I_end,
                       opts.of.interest$R_end,
                       opts.of.interest$D_end,
                       opts.of.interest$T_end,
                       opts.of.interest$TR_end,
                       opts.of.interest$TD_end)
    names(opts) <- c("R0","E_0","I_0","R_0","D_0","T_0","TR_0","TD_0")
    new.ops <- SEIR.against.real.data.multiple(simulation.start,simulation.end,opts,data,mod,opts.of.interest$reported.end)
    for(j in 1:length(new.ops)){
      save.opts[[length(save.opts)+1]] <- last.step[[i]]
      save.opts[[length(save.opts)]] <- rbind(save.opts[[length(save.opts)]],new.ops[[j]])
    }
  }
  return(save.opts)
}

modelSEIR <- function(data,opts.R0,opts.E_0,opts.I_0){
  
  names(data)[which(names(data)!="Date" & names(data)!="Underreporting")] <- "cases.observed"
  data$Date <- as.character(data$Date)
  data$Date <- gsub("\\.","-",as.character(data$Date))
  data$Date <- as.Date(data$Date,format="%d-%m-%Y")
  data$day <- as.numeric(data$Date-min(data$Date))+1
  data$days.to.today <- as.numeric(as.Date(format(Sys.time(), "%Y-%m-%d"))-data$Date)
  
  f <- data
  names(f)[which(names(f)=="day")] <- "time"
  f$Underreporting <- logit(1/f$Underreporting )
  mod <- lm(Underreporting~time^2,data=f)
  
  data <- data[!is.na(data$cases.observed),]
  
  today <- as.Date(format(Sys.time(), "%Y-%m-%d"))
  days.to.today <- min(data$days.to.today)
  days <- data$day[order(data$day)]
  
  start.date <- max(days)-floor(max(days)/days.break)*days.break+1
  while(sum(days<start.date+days.break-1)==0){
    start.date <- start.date+days.break
  }
  
  
  opts <- expand.grid(opts.R0,opts.E_0,opts.I_0)
  names(opts) <- c("R0","E_0","I_0")
  opts$R_0 <- data$cases.observed[data$day==min(data$day[data$day>=start.date])]-opts$E_0
  opts$D_0 <- 0
  opts$T_0 <- 0
  opts$TR_0 <- data$cases.observed[data$day==min(data$day[data$day>=start.date])]-opts$E_0
  opts$TD_0 <- 0
  opts <- opts[opts$R_0>=0,]
  
  step1 <- SEIR.against.real.data.multiple(start.date,start.date+days.break-1,opts,data,mod,0)
  step1 <- reduce.options(step1)
  
  first.day <- start.date+days.break
  while(first.day<max(days)){
    print("--")
    print(first.day)
    print(max(days))
    print(length(step1))
    step1 <- SEIR.second.step.against.real.data.multiple(first.day-1,first.day+days.break-1,step1,opts.R0,data,mod)
    step1 <- reduce.options(step1)
    first.day <- first.day+days.break
  }
  
  params <- step1
  for(i in 1:length(params)) params[[i]]$index <- i
  
  collapsed.final.opts <- calc.prob(params)
  collapsed.final.opts <- collapsed.final.opts[,dif:=NULL]
  
  params <- rbindlist(params)
  params <- merge(params,collapsed.final.opts,by=c("index"))
  params$est.end <- params$I_end+params$R_end+params$D_end
  
  results <- "a"
  old <- current <- params[1,]
  
  print(nrow(params))
  for(i in 1:(nrow(params)+1)){
    if(i/nrow(params)*100==floor(i/nrow(params)*100/1)*1){
      cat(i/nrow(params)*100)
      cat("\n")
    }
    old <- current
    if(i<=nrow(params))current <- params[i,]
    reported.end <- old$reported.end
    if(i==1) reported.end<-0
    
    new.run <- FALSE
    if(old$index!=current$index | i==(nrow(params)+1)){
      sirmodel <- SEIR(old$d.end,old$d.end+days.to.today+6*28,
                       R0=old$R0,
                       time.incubation=time.incubation,
                       time.infectedL=time.infectedL,
                       time.infectedD=time.infectedD,
                       cfr=cfr,
                       E_0=old$E_end,
                       I_0=old$I_end,
                       R_0=old$R_end,
                       D_0=old$D_end,
                       T_0=old$T_end,
                       TR_0=old$TR_end,
                       TD_0=old$TD_end,
                       mod,reported.end)
      sirmodel <- sirmodel[-1,]
      sirmodel$p <- old$p
      sirmodel$index <- old$index
      sirmodel$R0 <- old$R0
      sirmodel$type <- "projection"
      
      results <- rbind(results,sirmodel)
      
      sirmodel <- SEIR(old$d.end,old$d.end+days.to.today+6*28,
                       R0=old$R0,
                       time.incubation=time.incubation,
                       time.infectedL=time.infectedL,
                       time.infectedD=time.infectedD,
                       cfr=cfr,
                       E_0=old$E_end,
                       I_0=old$I_end,
                       R_0=old$R_end,
                       D_0=old$D_end,
                       T_0=old$T_end,
                       TR_0=old$TR_end,
                       TD_0=old$TD_end,
                       mod,reported.end,
                       referenceday=old$d.end+days.to.today,
                       reduced.day=old$d.end+days.to.today,
                       reduced.scenario=TRUE)
      sirmodel <- sirmodel[-1,]
      sirmodel$p <- old$p
      sirmodel$index <- old$index
      sirmodel$R0 <- old$R0
      sirmodel$type <- "cdcreduction_today_5dayreport"
      
      results <- rbind(results,sirmodel)
      
      sirmodel <- SEIR(old$d.end,old$d.end+days.to.today+6*28,
                       R0=old$R0,
                       time.incubation=time.incubation,
                       time.infectedL=time.infectedL,
                       time.infectedD=time.infectedD,
                       cfr=cfr,
                       E_0=old$E_end,
                       I_0=old$I_end,
                       R_0=old$R_end,
                       D_0=old$D_end,
                       T_0=old$T_end,
                       TR_0=old$TR_end,
                       TD_0=old$TD_end,
                       mod,reported.end,
                       referenceday=old$d.end+days.to.today,
                       reduced.day=old$d.end+days.to.today+28,
                       reduced.scenario=TRUE)
      sirmodel <- sirmodel[-1,]
      sirmodel$p <- old$p
      sirmodel$index <- old$index
      sirmodel$R0 <- old$R0
      sirmodel$type <- "cdcreduction_5dayreport_28"
      
      results <- rbind(results,sirmodel)
      
      sirmodel <- SEIR(old$d.end,old$d.end+days.to.today+6*28,
                       R0=old$R0,
                       time.incubation=time.incubation,
                       time.infectedL=time.infectedL,
                       time.infectedD=time.infectedD,
                       cfr=cfr,
                       E_0=old$E_end,
                       I_0=old$I_end,
                       R_0=old$R_end,
                       D_0=old$D_end,
                       T_0=old$T_end,
                       TR_0=old$TR_end,
                       TD_0=old$TD_end,
                       mod,reported.end,
                       referenceday=old$d.end+days.to.today,
                       reduced.day=old$d.end+days.to.today+2*28,
                       reduced.scenario=TRUE)
      sirmodel <- sirmodel[-1,]
      sirmodel$p <- old$p
      sirmodel$index <- old$index
      sirmodel$R0 <- old$R0
      sirmodel$type <- "cdcreduction_5dayreport_56"
      
      results <- rbind(results,sirmodel)
      
      #
      
      sirmodel <- SEIR(old$d.end,old$d.end+days.to.today+6*28,
                       R0=old$R0,
                       time.incubation=time.incubation,
                       time.infectedL=time.infectedL,
                       time.infectedD=time.infectedD,
                       cfr=cfr,
                       E_0=old$E_end,
                       I_0=old$I_end,
                       R_0=old$R_end,
                       D_0=old$D_end,
                       T_0=old$T_end,
                       TR_0=old$TR_end,
                       TD_0=old$TD_end,
                       mod,reported.end,
                       referenceday=old$d.end+days.to.today,
                       reduced.day=old$d.end+days.to.today,
                       reduced.scenario=TRUE,
                       time.to.report=3
                       )
      sirmodel <- sirmodel[-1,]
      sirmodel$p <- old$p
      sirmodel$index <- old$index
      sirmodel$R0 <- old$R0
      sirmodel$type <- "cdcreduction_today_3dayreport"
      
      results <- rbind(results,sirmodel)
      
      sirmodel <- SEIR(old$d.end,old$d.end+days.to.today+6*28,
                       R0=old$R0,
                       time.incubation=time.incubation,
                       time.infectedL=time.infectedL,
                       time.infectedD=time.infectedD,
                       cfr=cfr,
                       E_0=old$E_end,
                       I_0=old$I_end,
                       R_0=old$R_end,
                       D_0=old$D_end,
                       T_0=old$T_end,
                       TR_0=old$TR_end,
                       TD_0=old$TD_end,
                       mod,reported.end,
                       referenceday=old$d.end+days.to.today,
                       reduced.day=old$d.end+days.to.today+28,
                       reduced.scenario=TRUE,
                       time.to.report=3
                       )
      sirmodel <- sirmodel[-1,]
      sirmodel$p <- old$p
      sirmodel$index <- old$index
      sirmodel$R0 <- old$R0
      sirmodel$type <- "cdcreduction_3dayreport_28"
      
      results <- rbind(results,sirmodel)
      
      sirmodel <- SEIR(old$d.end,old$d.end+days.to.today+6*28,
                       R0=old$R0,
                       time.incubation=time.incubation,
                       time.infectedL=time.infectedL,
                       time.infectedD=time.infectedD,
                       cfr=cfr,
                       E_0=old$E_end,
                       I_0=old$I_end,
                       R_0=old$R_end,
                       D_0=old$D_end,
                       T_0=old$T_end,
                       TR_0=old$TR_end,
                       TD_0=old$TD_end,
                       mod,reported.end,
                       referenceday=old$d.end+days.to.today,
                       reduced.day=old$d.end+days.to.today+2*28,
                       reduced.scenario=TRUE,
                       time.to.report=3
                       )
      sirmodel <- sirmodel[-1,]
      sirmodel$p <- old$p
      sirmodel$index <- old$index
      sirmodel$R0 <- old$R0
      sirmodel$type <- "cdcreduction_3dayreport_56"
      
      results <- rbind(results,sirmodel)
      
      reported.end <- 0
      new.run <- TRUE
    } 
    
    if(i!=(nrow(params)+1)){
    
      sirmodel <- SEIR(current$d.start,current$d.end,
                       R0=current$R0,
                       time.incubation=time.incubation,
                       time.infectedL=time.infectedL,
                       time.infectedD=time.infectedD,
                       cfr=cfr,
                       E_0=current$E_0,
                       I_0=current$I_0,
                       R_0=current$R_0,
                       D_0=current$D_0,
                       T_0=current$T_0,
                       TR_0=current$TR_0,
                       TD_0=current$TD_0,
                       mod,reported.end)
    
      if(!new.run & i!=1) sirmodel <- sirmodel[-1,]
      sirmodel$p <- current$p
      sirmodel$index <- current$index
      sirmodel$R0 <- current$R0
      sirmodel$type <- "historic"
      if(is.character(results)){
        results <- sirmodel
      } else results <- rbind(results,sirmodel)
    }
  }
  ############################################
  results <- data.table(results)
  results.p <- results[,list(p=sum(p)),by=list(type,time)]
  setnames(results.p,c("type","time","pdenominator"))
  adj.results <- merge(results,results.p,by=c("type","time"))
  adj.results$adjp <- adj.results$p/adj.results$pdenominator
  
  collapsed <- adj.results[,
                      list(R0=sum(R0*adjp),
                           R02=sum(R0^2*adjp),
                           cases.estimated=sum(cases.estimated*adjp),
                           cases.estimated2=sum(cases.estimated^2*adjp),
                           cases.reported=sum(cases.reported*adjp),
                           cases.reported2=sum(cases.reported^2*adjp),
  				 E=sum(E*adjp),
  				 E2=sum(E^2*adjp),
  				 I=sum(I*adjp),
  				 I2=sum(I^2*adjp),
           T=sum(T*adjp),
           T2=sum(T^2*adjp)
  				 ),
                      by=list(time,type)]
  
  collapsed$R0.sd <- sqrt(collapsed$R02-collapsed$R0^2)
  collapsed$R0.sd[is.nan(collapsed$R0.sd)] <- 0
  
  collapsed$cases.estimated.sd <- sqrt(collapsed$cases.estimated2-collapsed$cases.estimated^2)
  collapsed$cases.estimated.sd[is.nan(collapsed$cases.estimated.sd)] <- 0

  collapsed$cases.reported.sd <- sqrt(collapsed$cases.reported2-collapsed$cases.reported^2)
  collapsed$cases.reported.sd[is.nan(collapsed$cases.reported.sd)] <- 0

  collapsed$E.sd <- sqrt(collapsed$E2-collapsed$E^2)
  collapsed$E.sd[is.nan(collapsed$E.sd)] <- 0

  collapsed$I.sd <- sqrt(collapsed$I2-collapsed$I^2)
  collapsed$I.sd[is.nan(collapsed$I.sd)] <- 0

  collapsed$T.sd <- sqrt(collapsed$T2-collapsed$T^2)
  collapsed$T.sd[is.nan(collapsed$T.sd)] <- 0

  collapsed <- collapsed[,R02:=NULL]
  collapsed <- collapsed[,cases.estimated2:=NULL]
  collapsed <- collapsed[,E2:=NULL]
  collapsed <- collapsed[,I2:=NULL]
  collapsed <- collapsed[,T2:=NULL]
  setnames(collapsed,1,"day")
  
  collapsed$cases.estimated.l95 <- collapsed$cases.estimated-1.96*collapsed$cases.estimated.sd
  collapsed$cases.estimated.u95 <- collapsed$cases.estimated+1.96*collapsed$cases.estimated.sd
  
  collapsed$cases.reported.l95 <- collapsed$cases.reported-1.96*collapsed$cases.reported.sd
  collapsed$cases.reported.u95 <- collapsed$cases.reported+1.96*collapsed$cases.reported.sd
  
  collapsed$R0.l95 <- collapsed$R0-1.96*collapsed$R0.sd
  collapsed$R0.u95 <- collapsed$R0+1.96*collapsed$R0.sd
  
  collapsed$E.l95 <- collapsed$E-1.96*collapsed$E.sd
  collapsed$E.u95 <- collapsed$E+1.96*collapsed$E.sd
  
  collapsed$I.l95 <- collapsed$I-1.96*collapsed$I.sd
  collapsed$I.u95 <- collapsed$I+1.96*collapsed$I.sd
  
  collapsed$T.l95 <- collapsed$T-1.96*collapsed$T.sd
  collapsed$T.u95 <- collapsed$T+1.96*collapsed$T.sd
  
  collapsed$days.to.today <- collapsed$day-max(days)-days.to.today
  collapsed$date <- today+collapsed$days.to.today
  
  collapsed <- merge(collapsed,data[,c("day","cases.observed")],by="day",all=TRUE)
  
	f <- data.frame(time=collapsed$day)
	collapsed$reporting.quotient <- inv.logit(predict(mod,f))
  collapsed$reporting.quotient[collapsed$days.to.today>=0 & collapsed$type=="cdcreduction_today"] <- 0.7
  collapsed$reporting.quotient[collapsed$days.to.today>=28 & collapsed$type=="cdcreduction_28"] <- 0.7
  collapsed$reporting.quotient[collapsed$days.to.today>=2*28 & collapsed$type=="cdcreduction_56"] <- 0.7

  return(collapsed)
}

data <- read.csv("data.csv",header=TRUE,sep=";")

data.models <- list()
data.models[["liberia"]] <- data[,c("Underreporting","Date","LiberiaCases")]
data.models[["sierra"]] <- data[,c("Underreporting","Date","SierraCases")]
data.models[["guinea"]] <- data[,c("Underreporting","Date","GuineaCases")]

opts.R0 <- seq(1.2,2.4,0.1)
opts.E_0 <- seq(2,82,10)#c(seq(2,20,4),seq(30,200,20))
opts.I_0 <- seq(2,82,10)#c(seq(2,20,4),seq(30,200,20))

data <- data.models[["sierra"]]

results.models <- list()
for(i in names(data.models)){
	print(i)
	results.models[[i]] <- modelSEIR(data.models[[i]],opts.R0,opts.E_0,opts.I_0)
}

try(results.models[["total"]]<-NULL,TRUE)
aggregated <- rbindlist(results.models)
aggregated <- aggregated[!is.na(aggregated$days.to.today),]
aggregated <- aggregated[,
	list(
	  reporting.quotient=mean(reporting.quotient),
	  cases.estimated=sum(cases.estimated),
	  cases.estimated.sd=sum(cases.estimated.sd^2),
	  cases.reported=sum(cases.reported),
	  cases.reported.sd=sum(cases.reported.sd^2),
		cases.observed=sum(cases.observed),
				 E=sum(E),
				 E.sd=sum(E.sd^2),
				 I=sum(I),
				 I.sd=sum(I.sd^2),
         T=sum(T),
         T.sd=sum(T.sd^2)),
	by=list(days.to.today,type)]



aggregated$cases.estimated.sd <- sqrt(aggregated$cases.estimated.sd)
aggregated$cases.estimated.l95 <- aggregated$cases.estimated-1.96*aggregated$cases.estimated.sd
aggregated$cases.estimated.u95 <- aggregated$cases.estimated+1.96*aggregated$cases.estimated.sd

aggregated$cases.reported.sd <- sqrt(aggregated$cases.reported.sd)
aggregated$cases.reported.l95 <- aggregated$cases.reported-1.96*aggregated$cases.reported.sd
aggregated$cases.reported.u95 <- aggregated$cases.reported+1.96*aggregated$cases.reported.sd

aggregated$E.sd <- sqrt(aggregated$E.sd)
aggregated$E.l95 <- aggregated$E-1.96*aggregated$E.sd
aggregated$E.u95 <- aggregated$E+1.96*aggregated$E.sd

aggregated$I.sd <- sqrt(aggregated$I.sd)
aggregated$I.l95 <- aggregated$I-1.96*aggregated$I.sd
aggregated$I.u95 <- aggregated$I+1.96*aggregated$I.sd

aggregated$T.sd <- sqrt(aggregated$T.sd)
aggregated$T.l95 <- aggregated$T-1.96*aggregated$T.sd
aggregated$T.u95 <- aggregated$T+1.96*aggregated$T.sd

results.models[["total"]] <- aggregated

for(i in names(results.models)){
  results.models[[i]] <- results.models[[i]][-(1:28),]
}

for(i in names(results.models)){
  results.models[[i]]$cdc25p.ETU <- 0.25*(results.models[[i]]$I+results.models[[i]]$T)
  results.models[[i]]$cdc25p.ETU.sd <- 0.25*sqrt(results.models[[i]]$I.sd^2+results.models[[i]]$T.sd^2)
  results.models[[i]]$cdc25p.ETU.l95 <- results.models[[i]]$cdc25p.ETU-1.96*results.models[[i]]$cdc25p.ETU.sd
  results.models[[i]]$cdc25p.ETU.u95 <- results.models[[i]]$cdc25p.ETU+1.96*results.models[[i]]$cdc25p.ETU.sd

  results.models[[i]]$cdc45p.reduced <- 0.45*(results.models[[i]]$I+results.models[[i]]$T)
  results.models[[i]]$cdc45p.reduced.sd <- 0.45*sqrt(results.models[[i]]$I.sd^2+results.models[[i]]$T.sd^2)
  results.models[[i]]$cdc45p.reduced.l95 <- results.models[[i]]$cdc45p.reduced-1.96*results.models[[i]]$cdc45p.reduced.sd
  results.models[[i]]$cdc45p.reduced.u95 <- results.models[[i]]$cdc45p.reduced+1.96*results.models[[i]]$cdc45p.reduced.sd
}

plotting.models <- results.models
for(j in names(plotting.models)){
  plotting.models[[j]]$label <- ""
  plotting.models[[j]]$label[plotting.models[[j]]$type=="historic"] <- "Historical"
  plotting.models[[j]]$label[plotting.models[[j]]$type=="projection"] <- "Projection"
  plotting.models[[j]]$label[plotting.models[[j]]$type=="cdcreduction_5dayreport_56"] <- "Target in 56 days\n5 days to hospital"
  plotting.models[[j]]$label[plotting.models[[j]]$type=="cdcreduction_5dayreport_28"] <- "Target in 28 days\n5 days to hospital"
  plotting.models[[j]]$label[plotting.models[[j]]$type=="cdcreduction_today_5dayreport"] <- "Target in 0 days\n5 days to hospital"
  plotting.models[[j]]$label[plotting.models[[j]]$type=="cdcreduction_3dayreport_56"] <- "Target in 56 days\n3 days to hospital"
  plotting.models[[j]]$label[plotting.models[[j]]$type=="cdcreduction_3dayreport_28"] <- "Target in 28 days\n3 days to hospital"
  plotting.models[[j]]$label[plotting.models[[j]]$type=="cdcreduction_today_3dayreport"] <- "Target in 0 days\n3 days to hospital"
  plotting.models[[j]]$label <- factor(plotting.models[[j]]$label,levels=c(
    "Historical",
    "Projection",
    "Target in 56 days\n5 days to hospital",
    "Target in 28 days\n5 days to hospital",
    "Target in 0 days\n5 days to hospital",
    "Target in 56 days\n3 days to hospital",
    "Target in 28 days\n3 days to hospital",
    "Target in 0 days\n3 days to hospital"
    ))
}

labeled.data <- plotting.models

ybreaks <- c(1000,2000,5000,10000,20000,40000,100000,200000,400000,1000000,2000000,4000000,10000000,20000000,40000000)
ylabs <- c(paste(ybreaks[1:9]/1000,"K",sep=""),paste(ybreaks[10:length(ybreaks)]/1000000,"M",sep=""))

plotting.models <- labeled.data
for(j in names(plotting.models)){
  plotting.models[[j]]$cases.estimated.l95[plotting.models[[j]]$days.to.today>58] <- NA
  plotting.models[[j]]$cases.estimated.l95[plotting.models[[j]]$days.to.today>58] <- NA
  xbreaks <- seq(-28*5,28*6,28)
  #text <- collapsed[nrow(collapsed),]
  q <- ggplot(plotting.models[[j]][plotting.models[[j]]$label!="Historical",],aes(x=days.to.today))
  for(i in ybreaks) q <- q + geom_hline(yintercept=log2(i),lty=3,lwd=2)
  for(i in xbreaks) q <- q + geom_vline(xintercept=i,lty=3,lwd=2)
  q <- q + geom_ribbon(data=plotting.models[[j]][plotting.models[[j]]$label=="Projection",],mapping=aes(ymin=log2(cases.estimated.l95),ymax=log2(cases.estimated.u95),fill=label),
                       alpha=0.4)
  q <- q + geom_line(aes(y=log2(cases.estimated),col=label),lwd=3)
  
  q <- q + scale_colour_brewer("",palette="Set1")
  q <- q + scale_fill_brewer("",palette="Set1")
  q <- q + guides(colour = guide_legend(override.aes = list(size=30),reverse=FALSE))
  q <- q + guides(fill=FALSE)
  
  q <- q + labs(title=paste(proper(j)," estimates of cumulative cases corrected for underreporting as of ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),"\n",sep=""))
  q <- q + scale_x_continuous(paste("\nDays from ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),sep=""),breaks=xbreaks,lim=c(0,168))
  q <- format_plot(q,6,10)
  q <- q + scale_y_continuous("Cumulative cases\n",breaks=log2(ybreaks),labels=ylabs )
  png(paste(today.character,"/",today.character,"_intervention_log_estimates_",j,".png",sep=""),width=3508,height=2480)
  print(q)
  dev.off()
}


plotting.models <- labeled.data
ybreaks <- c(1000,2000,5000,10000,20000,40000,100000,200000)
ylabs <- paste(ybreaks/1000,"K",sep="")

for(j in names(plotting.models)){
  plotting.models[[j]] <- plotting.models[[j]][plotting.models[[j]]$label=="Historical",]
  xbreaks <- seq(-28*5,0,28)
  
  q <- ggplot(plotting.models[[j]],aes(x=days.to.today))
  for(i in ybreaks) q <- q + geom_hline(yintercept=log2(i),lty=3,lwd=2)
  for(i in xbreaks) q <- q + geom_vline(xintercept=i,lty=3,lwd=2)
  q <- q + geom_point(aes(y=log2(cases.observed)),size=20)
  
  q <- q + geom_ribbon(aes(ymin=log2(cases.reported.l95),ymax=log2(cases.reported.u95),fill="blue"),
                       alpha=0.4)
  q <- q + geom_line(aes(y=log2(cases.reported),col="blue"),lwd=3)
  
  q <- q + geom_ribbon(aes(ymin=log2(cases.estimated.l95),ymax=log2(cases.estimated.u95),fill="red"),
                       alpha=0.4)
  q <- q + geom_line(aes(y=log2(cases.estimated),col="red"),lwd=3)
  
  q <- q + scale_colour_manual("",values=c("blue"="blue","red"="red"),labels=c("Reported","Estimated"))
  q <- q + scale_fill_manual("",values=c("blue"="blue","red"="red"),labels=c("Reported","Estimated"))
  
  q <- q + guides(colour = guide_legend(override.aes = list(size=30),reverse=TRUE))
  q <- q + guides(fill=FALSE)
  
  q <- q + labs(title=paste(proper(j)," estimates of cumulative cases as of ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),"\n",sep=""))
  q <- q + scale_x_continuous(paste("\nDays from ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),sep=""),breaks=xbreaks)
  q <- format_plot(q,6,10)
  q <- q + scale_y_continuous("Cumulative cases\n",breaks=log2(ybreaks),labels=ylabs )
  png(paste(today.character,"/",today.character,"_log_validation_",j,".png",sep=""),width=3508,height=2480)
  print(q)
  dev.off()
}

ybreaks <- c(5,10,20,50,100,200,500,1000,2000,5000,10000,20000,40000,100000,200000)
ylabs <- c(ybreaks[1:7],paste(ybreaks[8:length(ybreaks)]/1000,"K",sep=""))
xbreaks <- seq(-28*5,6*28,28)
new.day <- labeled.data

store.data <- list()
for(j in names(new.day)){
  x <- new.day[[j]]
  l <- list()
  for(i in unique(x$label)){
    l[[i]] <- x[x$label==i,]
  }
  for(i in unique(x$label)){
    if(i!="Historical") l[[i]] <- rbind(as.data.frame(l[["Historical"]])[,names(l[[i]])],l[[i]])
    l[[i]]$new.cases.estimated <- l[[i]]$cases.estimated
    l[[i]]$new.cases.estimated[2:nrow(l[[i]])] <- l[[i]]$cases.estimated[2:nrow(l[[i]])]-l[[i]]$cases.estimated[1:(nrow(l[[i]])-1)]
    l[[i]]$new.cases.estimated[1] <- NA
    
    l[[i]]$new.cases.reported <- l[[i]]$cases.reported
    l[[i]]$new.cases.reported[2:nrow(l[[i]])] <- l[[i]]$cases.reported[2:nrow(l[[i]])]-l[[i]]$cases.reported[1:(nrow(l[[i]])-1)]
    l[[i]]$new.cases.reported[1] <- NA
    if(i!="Historical") l[[i]] <- l[[i]][l[[i]]$label!="Historical",]
  }
  
  x <- l[["Historical"]]
  for(i in names(l)){
    if(i=="Historical") next
    x <- rbind(x,l[[i]])
  }
  store.data[[j]] <- x
  
  q <- ggplot(x,aes(x=days.to.today))
  for(i in ybreaks) q <- q + geom_hline(yintercept=log2(i),lty=3,lwd=2)
  for(i in xbreaks) q <- q + geom_vline(xintercept=i,lty=3,lwd=2)
#  q <- q + geom_line(aes(y=log2(new.cases.reported),col="blue"),lwd=3)
  q <- q + geom_line(aes(y=log2(new.cases.estimated),col=label),lwd=3)
  
  q <- q + scale_colour_brewer("",palette="Set1")
  q <- q + scale_fill_brewer("",palette="Set1")
  q <- q + guides(colour = guide_legend(override.aes = list(size=30),reverse=FALSE))
  q <- q + guides(fill=FALSE)
 # q <- q + scale_colour_manual("",values=c("blue"="blue","red"="red"),labels=c("Reported","Estimated"))
  #q <- q + guides(colour = guide_legend(override.aes = list(size=30),reverse=TRUE))
  
  q <- q + labs(title=paste(proper(j)," estimates of new daily cases corrected for underreporting as of ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),"\n",sep=""))
  q <- q + scale_x_continuous(paste("\nDays from ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),sep=""),breaks=xbreaks)
  q <- format_plot(q,6,10)
  q <- q + scale_y_continuous("New daily cases\n",breaks=log2(ybreaks),labels=ylabs)
  png(paste(today.character,"/",today.character,"_new_cases_",j,".png",sep=""),width=3508,height=2480)
  print(q)
  dev.off()
}


###################################



new.dayx <- labeled.data
for(j in names(results.models)){
ybreaks <- c(100,200,500,1000,2000,5000,10000,20000,40000,100000)
ylabs <- paste(ybreaks/1000,"K",sep="")
xbreaks <- seq(-28*5,28*2,28)
#text <- collapsed[nrow(collapsed),]

new.dayx[[j]]$cdc25p.ETU.l95[new.dayx[[j]]$cdc25p.ETU.l95<100] <- 100
new.dayx[[j]]$cdc45p.reduced.l95[new.dayx[[j]]$cdc45p.reduced.l95<100] <- 100

q <- ggplot(new.dayx[[j]][new.dayx[[j]]$label=="Historical"|new.dayx[[j]]$label=="Projection",],aes(x=days.to.today))
for(i in ybreaks) q <- q + geom_hline(yintercept=log2(i),lty=3,lwd=2)
for(i in xbreaks) q <- q + geom_vline(xintercept=i,lty=3,lwd=2)

q <- q + geom_ribbon(aes(ymin=log2(cdc25p.ETU.l95),ymax=log2(cdc25p.ETU.u95),fill=label),
                     alpha=0.4)
q <- q + geom_line(aes(y=log2(cdc25p.ETU),col=label,lty="ETU"),lwd=3)

q <- q + geom_ribbon(aes(ymin=log2(cdc45p.reduced.l95),ymax=log2(cdc45p.reduced.u95),fill=label),
                     alpha=0.4)
#q <- q + geom_ribbon(aes(ymin=log2(cdc45p.reduced.l95),ymax=log2(cdc45p.reduced.u95),fill="red"),
#                    alpha=0.4)
q <- q + geom_line(aes(y=log2(cdc45p.reduced),col=label,lty="Reduced\ntransmission"),lwd=3)

q <- q + scale_linetype_manual("\nCDC targets",values=c("ETU"=1,"Reduced\ntransmission"=2),labels=c("ETU","Reduced\ntransmission"))
q <- q + scale_colour_brewer("",palette="Set1")
q <- q + scale_fill_brewer("",palette="Set1")
q <- q + guides(colour = guide_legend(override.aes = list(size=30),reverse=FALSE))
q <- q + guides(linetype = guide_legend(reverse=TRUE))
q <- q + guides(fill=FALSE)

q <- q + labs(title=paste(proper(j)," estimates of cases needing treatment to reach CDC 70% goals as of ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),"\n",sep=""))
q <- q + scale_x_continuous(paste("\nDays from ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),sep=""),breaks=xbreaks,lim=c(-56,56))
q <- format_plot(q,6,10)
q <- q + scale_y_continuous("",breaks=log2(ybreaks),labels=ylabs,lim=c(log2(100),log2(100000)))
png(paste(today.character,"/",today.character,"_CDC_limits_",j,".png",sep=""),width=3508,height=2480)
print(q)
dev.off()
}

for(j in names(new.day)){
  if(j!="total"){
	print("----------")
	print(j)
	print(new.day[[j]][1,])
  }
  to.print <- store.data[[j]]
  to.print <- to.print[to.print$label=="Historical" | to.print$label=="Projection",]
  to.print <- data.frame(to.print[order(-to.print$days.to.today),])
  to.print$date <- as.Date(today.character)+to.print$days.to.today
  to.print$cases.estimated <- format_est_ci(to.print$cases.estimated, to.print$cases.estimated.l95, to.print$cases.estimated.u95,0)
  to.print$cases.reported <- format_est_ci(to.print$cases.reported, to.print$cases.reported.l95, to.print$cases.reported.u95,0)
  to.print$new.cases.reported.treated <- round(to.print$new.cases.reported)
  to.print$new.cases.estimated <- round(to.print$new.cases.estimated)

  to.print$incubation.estimated <- format_est_ci(to.print$E, to.print$E.l95, to.print$E.u95,0)
  to.print$infectious.estimated <- format_est_ci(to.print$I, to.print$I.l95, to.print$I.u95,0)
  to.print$currently.treated.estimated <- format_est_ci(to.print$T, to.print$T.l95, to.print$T.u95,0)
  to.print$reporting.quotient <- paste(format_sig(to.print$reporting.quotient*100,0),"%",sep="")
  to.print$cdc25p.ETU <- format_est_ci(to.print$cdc25p.ETU ,to.print$cdc25p.ETU.l95,to.print$cdc25p.ETU.u95,0)
  to.print$cdc45p.reduced <- format_est_ci(to.print$cdc45p.reduced,to.print$cdc45p.reduced.l95,to.print$cdc45p.reduced.u95,0)

  to.print <- to.print[,c("date","days.to.today","cases.reported","new.cases.reported.treated","cases.estimated","new.cases.estimated","incubation.estimated","infectious.estimated","currently.treated.estimated",
	"cdc25p.ETU","cdc45p.reduced","reporting.quotient")]

  for(i in 1:ncol(to.print)) to.print[,i] <- as.character(to.print[,i])
  titles <- to.print[1,]
  titles$date <- "Date"
  titles$days.to.today <- "Days from today"
  titles$cases.reported <- "Cases (reported)"
  titles$new.cases.reported.treated <- "New daily cases (reported)"
  titles$cases.estimated <- "Cases (total)"
  titles$new.cases.estimated <- "New daily cases (total)"
  titles$incubation.estimated <- "Incubation period"
  titles$infectious.estimated <- "Infectious and unreported"
  titles$currently.treated.estimated <- "Undergoing treatment"
  titles$cdc25p.ETU <- "CDC 25% ETU target"
  titles$cdc45p.reduced <- "CDC 45% reduced transmission target"
  titles$reporting.quotient <- "Reporting quotient"
  if(j=="total"){
	to.save.table <- to.print[to.print$days.to.today==0 | to.print$days.to.today==28 | to.print$days.to.today==56,]
      to.save.table <- rbind(titles,to.save.table)
  }
  titles <- titles[,names(to.print)]
  to.print <- rbind(titles,to.print)
  write.table(na.omit(to.print),file=paste(today.character,"/",today.character,"_written_estimates_",j,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=";",qmethod="double")
}


for(j in names(new.day)){
  if(j!="total"){
    print("----------")
    print(j)
    print(new.day[[j]][1,])
  }
  to.print <- store.data[[j]]
  to.print <- to.print[to.print$label=="Historical" | to.print$label=="Projection",]
  to.print <- data.frame(to.print[order(-to.print$days.to.today),])
  to.print$date <- as.Date(today.character)+to.print$days.to.today
  to.print$need.treatment <- round(to.print$cdc25p.ETU + to.print$cdc45p.reduced)
  to.print$needed.administrator.coordination <- round(to.print$need.treatment/20000*22)
  to.print$needed.doctor.etu <- round(to.print$need.treatment/20000*45)
  to.print$needed.doctor.coordination <- round(to.print$need.treatment/20000*11)
  to.print$needed.epidemiologist.contact.tracing <- round(to.print$need.treatment/20000*25)
  to.print$needed.epidemiologist.coordination <- round(to.print$need.treatment/20000*33)
  to.print$needed.labscientist.labs <- round(to.print$need.treatment/20000*67)
  to.print$needed.logistician.burials <- round(to.print$need.treatment/20000*33)
  to.print$needed.logistician.etu <- round(to.print$need.treatment/20000*89)
  to.print$needed.logistician.coordination <- round(to.print$need.treatment/20000*33)
  to.print$needed.nurse.etu <- round(to.print$need.treatment/20000*178)
  to.print$needed.phs.social <- round(to.print$need.treatment/20000*51)
  to.print$needed.phs.coordination <- round(to.print$need.treatment/20000*22)
  to.print$needed.datamanagement.coordination <- round(to.print$need.treatment/20000*28)
  to.print$total <- to.print$needed.administrator.coordination+
    to.print$needed.doctor.etu+
    to.print$needed.doctor.coordination+
    to.print$needed.epidemiologist.contact.tracing+
    to.print$needed.epidemiologist.coordination+
    to.print$needed.labscientist.labs+
    to.print$needed.logistician.burials+
    to.print$needed.logistician.etu+
    to.print$needed.logistician.coordination+
    to.print$needed.nurse.etu+
    to.print$needed.phs.social+
    to.print$needed.phs.coordination+
    to.print$needed.datamanagement.coordination
  
  to.print <- to.print[,c("date","days.to.today",
                          "need.treatment","needed.administrator.coordination","needed.doctor.etu","needed.doctor.coordination",
                          "needed.epidemiologist.contact.tracing","needed.epidemiologist.coordination","needed.labscientist.labs",
                          "needed.logistician.burials","needed.logistician.etu","needed.logistician.coordination",
                          "needed.nurse.etu","needed.phs.social","needed.phs.coordination","needed.datamanagement.coordination",
                          "total")]
  
  for(i in 1:ncol(to.print)) to.print[,i] <- as.character(to.print[,i])
  titles <- to.print[1,]
  titles$date <- "Date"
  titles$days.to.today <- "Days from today"
  titles$need.treatment <- "Caseload"
  titles$needed.administrator.coordination <- "Administrator (coordination)"
  titles$needed.doctor.etu <- "Doctor (ETU)"
  titles$needed.doctor.coordination <- "Doctor (coordination)"
  titles$needed.epidemiologist.contact.tracing <- "Epidemiologist (contact tracing)"
  titles$needed.epidemiologist.coordination <- "Epidemiologist (coordination)"
  titles$needed.labscientist.labs <- "Lab scientist (laboratories)"
  titles$needed.logistician.burials <- "Logistician (burials)"
  titles$needed.logistician.etu <- "Logistician (ETU)"
  titles$needed.logistician.coordination <- "Logistician (coordination)"
  titles$needed.nurse.etu <- "Nurse (ETU)"
  titles$needed.phs.social <- "Public health specialist (Social mobilization)"
  titles$needed.phs.coordination <- "Public health specialist (coordination)"
  titles$needed.datamanagement.coordination <- "Data management (coordination)"
  titles$total <- "Total personnel"
  
  if(j=="total"){
    to.save.table.persons <- to.print[to.print$days.to.today==0 | to.print$days.to.today==28 | to.print$days.to.today==56,]
    to.save.table.persons <- rbind(titles,to.save.table.persons)
  }
  titles <- titles[,names(to.print)]
  to.print <- rbind(titles,to.print)
  write.table(na.omit(to.print),file=paste(today.character,"/",today.character,"_persons_estimates_",j,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=";",qmethod="double")
}

for(j in names(new.day)){
  if(j!="total"){
    print("----------")
    print(j)
    print(new.day[[j]][1,])
  }
  to.print <- store.data[[j]]
  to.print <- to.print[to.print$label=="Historical" | to.print$label=="Projection",]
  to.print <- data.frame(to.print[order(-to.print$days.to.today),])
  to.print$date <- as.Date(today.character)+to.print$days.to.today
  to.print$cdc25p.ETU <- round(to.print$cdc25p.ETU)
  to.print$nurses.paramedics <- round(to.print$cdc25p.ETU/100*50)
  to.print$hygienists <- round(to.print$cdc25p.ETU/100*50)
  to.print$doctors <- round(to.print$cdc25p.ETU/100*10)
  to.print$support <- round(to.print$cdc25p.ETU)
  
  to.print$total <- to.print$nurses.paramedics+
    to.print$hygienists+
    to.print$doctors+
    to.print$support
  
  to.print <- to.print[,c("date","days.to.today",
                          "cdc25p.ETU","nurses.paramedics","hygienists","doctors",
                          "support","total")]
  
  for(i in 1:ncol(to.print)) to.print[,i] <- as.character(to.print[,i])
  titles <- to.print[1,]
  titles$date <- "Date"
  titles$days.to.today <- "Days from today"
  titles$cdc25p.ETU <- "ETU caseload"
  titles$nurses.paramedics <- "Nurses and paramedics"
  titles$hygienists <- "Hygienists and nurses aids"
  titles$doctors <- "Doctors"
  titles$support <- "Support staff"
  titles$total <- "Total"

  
  if(j=="total"){
    to.save.table.etu <- to.print[to.print$days.to.today==0 | to.print$days.to.today==28 | to.print$days.to.today==56,]
    to.save.table.etu <- rbind(titles,to.save.table.etu)
  }
  titles <- titles[,names(to.print)]
  to.print <- rbind(titles,to.print)
  write.table(na.omit(to.print),file=paste(today.character,"/",today.character,"_etu_estimates_",j,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=";",qmethod="double")
}

for(j in names(new.day)){
  if(j!="total"){
    print("----------")
    print(j)
    print(new.day[[j]][1,])
  }
  to.print <- store.data[[j]]
  to.print <- to.print[to.print$label=="Historical" | to.print$label=="Projection",]
  to.print <- data.frame(to.print[order(-to.print$days.to.today),])
  to.print$date <- as.Date(today.character)+to.print$days.to.today
  to.print$need.treatment <- round(to.print$cdc25p.ETU + to.print$cdc45p.reduced)
  
  to.print$needed.ppe.m3 <- round(to.print$need.treatment/20000*14893/(365.25/12),1)
  to.print$needed.ppe.tons <- round(to.print$need.treatment/20000*3095/(365.25/12),1)
  to.print$needed.bodybag.m3 <- round(to.print$need.treatment/20000*40/(365.25/12),1)
  to.print$needed.bodybag.tons <- round(to.print$need.treatment/20000*4/(365.25/12),1)
  to.print$needed.chlorine.m3 <- round(to.print$need.treatment/20000*16/(365.25/12),1)
  to.print$needed.chlorine.tons <- round(to.print$need.treatment/20000*16/(365.25/12),1)
  to.print$needed.medsup.m3 <- round(to.print$need.treatment/20000*230/(365.25/12),1)
  to.print$needed.medsup.tons <- round(to.print$need.treatment/20000*115/(365.25/12),1)
  to.print$needed.labsup.m3 <- round(to.print$need.treatment/20000*25/(365.25/12),1)
  to.print$needed.labsup.tons <- round(to.print$need.treatment/20000*12/(365.25/12),1)
  to.print$needed.total.m3 <- round(to.print$need.treatment/20000*15204/(365.25/12),1)
  to.print$needed.total.tons <- round(to.print$need.treatment/20000*3242/(365.25/12),1)
  
  
  to.print <- to.print[,c("date","days.to.today",
                          "needed.ppe.m3","needed.ppe.tons","needed.bodybag.m3","needed.bodybag.tons",
                          "needed.chlorine.m3","needed.chlorine.tons","needed.medsup.m3",
                          "needed.medsup.tons","needed.labsup.m3","needed.labsup.tons",
                          "needed.total.m3","needed.total.tons")]
  
  for(i in 1:ncol(to.print)) to.print[,i] <- as.character(to.print[,i])
  titles <- to.print[1,]
  titles$date <- "Date"
  titles$days.to.today <- "Days from today"
  titles$need.treatment <- "Caseload"
  titles$needed.ppe.m3 <- "PPE (m3)"
  titles$needed.ppe.tons <- "PPE (tons)"
  titles$needed.bodybag.m3 <- "Body bags (m3)"
  titles$needed.bodybag.tons <- "Body bags (tons)"
  titles$needed.chlorine.m3 <- "Chlorine (m3)"
  titles$needed.chlorine.tons <- "Chlorine (tons)"
  titles$needed.medsup.m3 <- "Medical supplies (m3)"
  titles$needed.medsup.tons <- "Medical supplies (tons)"
  titles$needed.labsup.m3 <- "Lab supplies (m3)"
  titles$needed.labsup.tons <- "Lab supplies (tons)"
  titles$needed.total.m3 <- "Total (m3)"
  titles$needed.total.tons <- "Total (tons)"
  titles <- titles[,names(to.print)]
  
  if(j=="total"){
    to.save.table.supplies <- to.print[to.print$days.to.today==0 | to.print$days.to.today==28 | to.print$days.to.today==56,]
    to.save.table.supplies <- rbind(titles,to.save.table.supplies)
  }
  titles <- titles[,names(to.print)]
  to.print <- rbind(titles,to.print)
  write.table(na.omit(to.print),file=paste(today.character,"/",today.character,"_supplies_estimates_",j,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=";",qmethod="double")
}


to.save <- store.data
for(j in names(store.data)){
	to.save[[j]]$date <- as.Date(today.character)+to.save[[j]]$days.to.today
}

saveRDS(to.save,file="raw_estimates.RDS")
saveRDS(to.save.table,file="table.RDS")
saveRDS(to.save.table.persons,file="table_persons.RDS")
saveRDS(to.save.table.etu,file="table_etu.RDS")
saveRDS(to.save.table.supplies,file="table_supplies.RDS")

rmarkdown::render("test.Rmd","word_document")
f <- paste(today.character,"/",today.character,"_ebola_model.docx",sep="")
try(file.remove(f),TRUE)
file.rename(from ="test.docx",  to =f)






