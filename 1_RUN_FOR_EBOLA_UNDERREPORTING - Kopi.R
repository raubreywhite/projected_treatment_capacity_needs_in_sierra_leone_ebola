##################################################################################
source("F:/_SM/SMAO/Felles/autostat/Libraries/r_tools/r_tools.R")
req(c("sfsmisc","deSolve","data.table"))

setwd("F:/_SM/_SM-Felles/Utbrudd/Ebola 2014/Ebola_models_updated")

days.break <- 28*2
time.incubation=11.4
time.infectedL=7.0
time.infectedD=7.0
cfr=0.6


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

SEIR.derivatives.function=function(t, x, vparameters){
  S = x[1]  # the value of S at time t
  E = x[2]
  I = x[3]  # the value of I at time t
  R = x[4]  # the value of R at time t
  D = x[5]

  if (E<0) E=0 # this is a cross check to ensure that we always have sensical values of E
  if (I<0) I=0 # this is a cross check to ensure that we always have sensical values of I
  
  with(as.list(vparameters),{
    npop = S+E+I+R   # the population size is always S+I+R because there are no births or deaths in the model
    dS = npop*mu - mu*S - beta*I/npop*S   # the derivative of S wrt time
    dE = beta*I/npop*S - (mu + a)*E  # the derivative of I wrt time
    dI = a*E - (gamma + mu)*I
    dR = (1-cfr)*gamma*I - mu*R                  # the derivative of R wrt time
    dD = cfr*gamma*I                  # the derivative of R wrt time
    out = c(dS,dE,dI,dR,dD)
    list(out)
  })
}

SEIR <- function(simulation.start,simulation.end,R0,time.incubation,time.infectedL,time.infectedD,cfr,
                 E_0=1,I_0=1,R_0=0,D_0=0,mod,reported.to.date){
  npop = 1000000
  S_0 = npop-E_0-I_0-R_0-D_0
  mu    = 0         # birth rate
  
  vt = seq(simulation.start,simulation.end,1)  # let's determine the values of S,I and R at times in vt
  a = 1/time.incubation         # approximate incubation period of influenza in days^{-1} 

  gammaL = 1/time.infectedL        # approximate average recovery period when living in days^{-1} 
  gammaD = 1/time.infectedD        # approximate average recovery period when dead in days^{-1} 
  gamma =  1/((1-cfr)/gammaL+cfr/gammaD)
  beta  = R0*(gamma+mu)*(a+mu)/a   # transmission rate of influenza is estimated by solving R0=beta/gamma for beta
  
  vparameters = c(gamma=gamma,gammaL=gammaL,gammaD=gammaD,beta=beta,mu=mu,cfr=cfr,a=a)
  inits = c(S=S_0,E=E_0,I=I_0,R=R_0,D=D_0)
  
  sirmodel = as.data.frame(lsoda(inits, vt, SEIR.derivatives.function, vparameters))
  
  sirmodel$Underreporting <- predict(mod,sirmodel)
  sirmodel$cases.estimated <- sirmodel$I+sirmodel$R+sirmodel$D
  sirmodel$cases.reported <- sirmodel$cases.estimated-sirmodel$cases.estimated[1]
  sirmodel$cases.reported[2:nrow(sirmodel)] <- sirmodel$cases.reported[2:nrow(sirmodel)]-sirmodel$cases.reported[1:(nrow(sirmodel)-1)]
  sirmodel$cases.reported <- sirmodel$cases.reported/sirmodel$Underreporting
  sirmodel$cases.reported <- reported.to.date+cumsum(sirmodel$cases.reported)
  sirmodel <- sirmodel[,-which(names(sirmodel)=="Underreporting")]
  
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
                       D_0=D_0,mod,reported.to.date)
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
  for(val in c(0.01,0.05,0.1,0.2,0.5)){
    collapsed.final.opts <- rbindlist(final.opts)
    collapsed.final.opts <- collapsed.final.opts[,
                                                 list(dif=mean(dif)),
                                                 by=list(index)]
    
    collapsed.final.opts$dif <- collapsed.final.opts$dif/(min(collapsed.final.opts$dif)*val)
    collapsed.final.opts$dif <- round(collapsed.final.opts$dif-min(collapsed.final.opts$dif),5)
    collapsed.final.opts$p <- exp(-collapsed.final.opts$dif)
    collapsed.final.opts$p <- round(collapsed.final.opts$p/sum(collapsed.final.opts$p),5)
    if(sum(collapsed.final.opts$p>0)>10) break
  }
  
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
                       opts.of.interest$D_end)
    names(opts) <- c("R0","E_0","I_0","R_0","D_0")
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
  
  for(i in 1:nrow(params)){
    old <- current
    current <- params[i,]
    reported.end <- old$reported.end
    if(i==1) reported.end<-0
    
    new.run <- FALSE
    if(old$index!=current$index){
      sirmodel <- SEIR(old$d.end,old$d.end+days.to.today+2*28,
                       R0=old$R0,
                       time.incubation=time.incubation,
                       time.infectedL=time.infectedL,
                       time.infectedD=time.infectedD,
                       cfr=cfr,
                       E_0=old$E_end,
                       I_0=old$I_end,
                       R_0=old$R_end,
                       D_0=old$D_end,mod,reported.end)
      sirmodel <- sirmodel[-1,]
      sirmodel$p <- old$p
      sirmodel$index <- old$index
      sirmodel$R0 <- old$R0
      
      results <- rbind(results,sirmodel)
      reported.end <- 0
      new.run <- TRUE
    } 
  
    sirmodel <- SEIR(current$d.start,current$d.end,
                     R0=current$R0,
                     time.incubation=time.incubation,
                     time.infectedL=time.infectedL,
                     time.infectedD=time.infectedD,
                     cfr=cfr,
                     E_0=current$E_0,
                     I_0=current$I_0,
                     R_0=current$R_0,
                     D_0=current$D_0,mod,reported.end)
  
    if(!new.run & i!=1) sirmodel <- sirmodel[-1,]
    sirmodel$p <- current$p
    sirmodel$index <- current$index
    sirmodel$R0 <- current$R0
    if(is.character(results)){
      results <- sirmodel
    } else results <- rbind(results,sirmodel)
  }
  
  results <- data.table(results)
  results.p <- results[,list(p=sum(p)),by=list(time)]
  setnames(results.p,c("time","pdenominator"))
  adj.results <- merge(results,results.p,by="time")
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
  				 I2=sum(I^2*adjp)
  				 ),
                      by=list(time)]
  
  collapsed$R0.sd <- sqrt(collapsed$R02-collapsed$R0^2)
  collapsed$cases.estimated.sd <- sqrt(collapsed$cases.estimated2-collapsed$cases.estimated^2)
  collapsed$cases.reported.sd <- sqrt(collapsed$cases.reported2-collapsed$cases.reported^2)
  collapsed$E.sd <- sqrt(collapsed$E2-collapsed$E^2)
  collapsed$I.sd <- sqrt(collapsed$I2-collapsed$I^2)
  collapsed <- collapsed[,R02:=NULL]
  collapsed <- collapsed[,cases.estimated2:=NULL]
  collapsed <- collapsed[,E2:=NULL]
  collapsed <- collapsed[,I2:=NULL]
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
  
  collapsed$days.to.today <- collapsed$day-max(days)-days.to.today
  collapsed$date <- today+collapsed$days.to.today
  
  collapsed <- merge(collapsed,data[,c("day","cases.observed")],by="day",all=TRUE)
  
  return(collapsed)
}

data <- read.csv("data.csv",header=TRUE,sep=";")

data.models <- list()
data.models[["liberia"]] <- data[,c("Underreporting","Date","LiberiaCases")]
data.models[["sierra"]] <- data[,c("Underreporting","Date","SierraCases")]
data.models[["guinea"]] <- data[,c("Underreporting","Date","GuineaCases")]

opts.R0 <- seq(0.8,2.4,0.1)
opts.E_0 <- seq(2,82,10)#c(seq(2,20,4),seq(30,200,20))
opts.I_0 <- seq(2,82,10)#c(seq(2,20,4),seq(30,200,20))

results.models <- list()
for(i in names(data.models)){
	results.models[[i]] <- modelSEIR(data.models[[i]],opts.R0,opts.E_0,opts.I_0)
}

try(results.models[["total"]]<-NULL,TRUE)
aggregated <- rbindlist(results.models)
aggregated <- aggregated[!is.na(aggregated$days.to.today),]
aggregated <- aggregated[,
	list(
	  cases.estimated=sum(cases.estimated),
	  cases.estimated.sd=sum(cases.estimated.sd^2),
	  cases.reported=sum(cases.reported),
	  cases.reported.sd=sum(cases.reported.sd^2),
		cases.observed=sum(cases.observed),
				 E=sum(E),
				 E.sd=sum(E.sd^2),
				 I=sum(I),
				 I.sd=sum(I.sd^2)),
	by=list(days.to.today)]



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

results.models[["total"]] <- aggregated

for(i in names(results.models)){
  results.models[[i]] <- results.models[[i]][-(1:28),]
}

#q <- ggplot(mapping=aes(x=day))
#q <- q + geom_line(data=collapsed,aes(y=R0),col="red")
#q

for(j in names(results.models)){
ybreaks <- c(50,100,1000,5000,10000,20000,40000,80000,160000)
xbreaks <- seq(-28*3,28*2,28)
#text <- collapsed[nrow(collapsed),]
q <- ggplot(results.models[[j]],aes(x=days.to.today))
for(i in ybreaks) q <- q + geom_hline(yintercept=log2(i),lty=2,lwd=2)
for(i in xbreaks) q <- q + geom_vline(xintercept=i,lty=2,lwd=2)
q <- q + geom_point(aes(y=log2(cases.observed)),size=10)

q <- q + geom_ribbon(aes(ymin=log2(cases.reported.l95),ymax=log2(cases.reported.u95)),
                     fill="blue",alpha=0.4)
q <- q + geom_line(aes(y=log2(cases.reported)),col="blue",lwd=3)

q <- q + geom_ribbon(aes(ymin=log2(cases.estimated.l95),ymax=log2(cases.estimated.u95)),
                    fill="red",alpha=0.4)
q <- q + geom_line(aes(y=log2(cases.estimated)),col="red",lwd=3)

q <- q + labs(title=paste(proper(j)," estimates as of ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),"\n",sep=""))
q <- q + scale_x_continuous(paste("\nDays from ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),sep=""),breaks=xbreaks)
q <- format_plot(q,6,6)
q <- q + scale_y_continuous("Cases\n",breaks=log2(ybreaks),labels=ybreaks)
png(paste(today.character,"/",today.character,"_log_estimates_",j,".png",sep=""),width=3508,height=2480)
print(q)
dev.off()
}

for(j in names(results.models)){
ybreaks <- c(50,1000,5000,10000,20000,40000,80000,160000)
xbreaks <- seq(-28*3,28*2,28)
#text <- collapsed[nrow(collapsed),]
q <- ggplot(results.models[[j]],aes(x=days.to.today))
for(i in ybreaks) q <- q + geom_hline(yintercept=i,lty=2,lwd=2)
for(i in xbreaks) q <- q + geom_vline(xintercept=i,lty=2,lwd=2)
q <- q + geom_ribbon(aes(ymin=cases.estimated.l95,ymax=cases.estimated.u95),
                    fill="red",alpha=0.4)
q <- q + geom_point(aes(y=cases.observed),size=10)
q <- q + geom_line(aes(y=cases.estimated),col="red",lwd=3)
q <- q + labs(title=paste(proper(j)," estimates as of ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),"\n",sep=""))
q <- q + scale_x_continuous(paste("\nDays from ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),sep=""),breaks=xbreaks)
q <- format_plot(q,6,6)
q <- q + scale_y_continuous("Cases\n",breaks=ybreaks,labels=ybreaks)
png(paste(today.character,"/",today.character,"_estimates_",j,".png",sep=""),width=3508,height=2480)
print(q)
dev.off()
}

for(j in names(results.models)){
if(j=="total") next
ybreaks <- seq(0.25,2,0.25)
xbreaks <- seq(-28*6,28*2,28)
q <- ggplot(results.models[[j]],aes(x=days.to.today))
for(i in ybreaks) q <- q + geom_hline(yintercept=i,lty=2,lwd=2)
q <- q + geom_ribbon(aes(ymin=R0.l95,ymax=R0.u95),
                    fill="red",alpha=0.4)
q <- q + geom_line(aes(y=R0),col="red",lwd=3)
q <- q + labs(title=paste(proper(j)," estimates as of ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),"\n",sep=""))
q <- q + scale_x_continuous(paste("\nDays from ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),sep=""),breaks=xbreaks)
q <- format_plot(q,6,6)
q <- q + scale_y_continuous("Cases\n",breaks=ybreaks,labels=ybreaks)
png(paste(today.character,"/",today.character,"_R0_estimates_",j,".png",sep=""),width=3508,height=2480)
print(q)
dev.off()
}

new.day <- results.models
for(j in names(new.day)){
  new.day[[j]]$new.cases.estimated <- new.day[[j]]$cases.estimated
  new.day[[j]]$new.cases.estimated[2:nrow(new.day[[j]])] <- new.day[[j]]$new.cases.estimated[2:nrow(new.day[[j]])]-new.day[[j]]$new.cases.estimated[1:(nrow(new.day[[j]])-1)]
  new.day[[j]]$new.cases.estimated[1] <- NA
  
  new.day[[j]]$new.cases.reported <- new.day[[j]]$cases.reported
  new.day[[j]]$new.cases.reported[2:nrow(new.day[[j]])] <- new.day[[j]]$new.cases.reported[2:nrow(new.day[[j]])]-new.day[[j]]$new.cases.reported[1:(nrow(new.day[[j]])-1)]
  new.day[[j]]$new.cases.reported[1] <- NA
  
  ybreaks <- c(5,10,25,50,100,200,500,1000,2000,5000,10000,20000,40000)
  xbreaks <- seq(-28*3,28*2,28)
  #text <- collapsed[nrow(collapsed),]
  q <- ggplot(new.day[[j]],aes(x=days.to.today))
  for(i in ybreaks) q <- q + geom_hline(yintercept=log2(i),lty=2,lwd=2)
  for(i in xbreaks) q <- q + geom_vline(xintercept=i,lty=2,lwd=2)
  q <- q + geom_line(aes(y=log2(new.cases.reported)),col="red",lwd=3)
  q <- q + geom_line(aes(y=log2(new.cases.estimated)),col="blue",lwd=3)
  q <- q + labs(title=paste(proper(j)," estimates as of ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),"\n",sep=""))
  q <- q + scale_x_continuous(paste("\nDays from ",as.character(as.Date(format(Sys.time(), "%Y-%m-%d"))),sep=""),breaks=xbreaks)
  q <- format_plot(q,6,6)
  q <- q + scale_y_continuous("New daily cases\n",breaks=log2(ybreaks),labels=ybreaks)
  png(paste(today.character,"/",today.character,"_new_cases_",j,".png",sep=""),width=3508,height=2480)
  print(q)
  dev.off()
}

for(j in names(new.day)){
  to.print <- new.day[[j]]
  to.print <- data.frame(to.print[order(-to.print$days.to.today),])
  to.print$date <- as.Date(today.character)+to.print$days.to.today
  to.print$cases.estimated <- format_est_ci(to.print$cases.estimated, to.print$cases.estimated.l95, to.print$cases.estimated.u95,0)
  to.print$cases.reported <- format_est_ci(to.print$cases.reported, to.print$cases.reported.l95, to.print$cases.reported.u95,0)
  to.print$new.cases.reported <- round(to.print$new.cases.reported)
  to.print$new.cases.estimated <- round(to.print$new.cases.estimated)

  to.print$incubation.estimated <- format_est_ci(to.print$E, to.print$E.l95, to.print$E.u95,0)
  to.print$infectious.estimated <- format_est_ci(to.print$I, to.print$I.l95, to.print$I.u95,0)


  to.print <- to.print[,c("date","days.to.today","cases.reported","new.cases.reported","cases.estimated","new.cases.estimated","incubation.estimated","infectious.estimated")]
  write.table(to.print,file=paste(today.character,"/",today.character,"_written_estimates_",j,".csv",sep=""),row.names=FALSE,sep=";",qmethod="double")
}


