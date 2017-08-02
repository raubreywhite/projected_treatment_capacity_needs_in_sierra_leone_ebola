##################################################################################
source("F:/_SM/SMAO/Felles/autostat/Libraries/r_tools/r_tools.R")
req(c("sfsmisc","deSolve","data.table"))

setwd("F:/_SM/_SM-Felles/Utbrudd/Ebola 2014/Ebola_models_updated")

days.break <- 28*1.5
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

SIRfunc_with_demographics=function(t, x, vparameters){
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
    dR = (1-cfr)*gammaL*I - mu*R                  # the derivative of R wrt time
    dD = cfr*gammaD*I                  # the derivative of R wrt time
    out = c(dS,dE,dI,dR,dD)
    list(out)
  })
}

SEIR <- function(simulation.start,simulation.end,R0,time.incubation,time.infectedL,time.infectedD,cfr,
                 E_0=1,I_0=1,R_0=0,D_0=0){
  npop = 1000000
  #E_0 = 1       # put one infected person in the population
  #I_0 = 0       
  S_0 = npop-E_0
  #R_0 = 0       # initially no one has recovered yet
  #D_0 = 0 # initially no one is dead
  
  mu    = 0         # birth rate
  
  vt = seq(simulation.start,simulation.end,1)  # let's determine the values of S,I and R at times in vt
  a = 1/time.incubation         # approximate incubation period of influenza in days^{-1} 

  gammaL = 1/time.infectedL        # approximate average recovery period when living in days^{-1} 
  gammaD = 1/time.infectedD        # approximate average recovery period when dead in days^{-1} 
  
  gamma = (1-cfr)*gammaL+cfr*gammaD
  
  
  beta  = R0*(gamma+mu)*(a+mu)/a   # transmission rate of influenza is estimated by solving R0=beta/gamma for beta
  
  
  vparameters = c(gamma=gamma,gammaL=gammaL,gammaD=gammaD,beta=beta,mu=mu,cfr=cfr,a=a)
  inits = c(S=S_0,E=E_0,I=I_0,R=R_0,D=D_0)
  
  sirmodel = as.data.frame(lsoda(inits, vt, SIRfunc_with_demographics, vparameters))
  sirmodel$cases.estimated <- sirmodel$I+sirmodel$R+sirmodel$D
  return(sirmodel)
}


fitSEIR <- function(simulation.start,simulation.end,opts,data){
  predictions <- list()
  opts$dif <- 0
  opts$obs.end <- 0
  for(i in 1:nrow(opts)){
    R0=opts[i,1]
    E_0=opts[i,2]
    I_0=opts[i,3]
    R_0=opts[i,4]
    D_0=opts[i,5]
    #max(data$days.to.today)+29
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
                     D_0=D_0)
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
    
    opts$dif[i] <- mean(((testing$cases.observed-testing$cases.estimated)^2))
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


trySEIR <- function(simulation.start,simulation.end,opts,data){
  retval <- fitSEIR(simulation.start,simulation.end,opts,data)
  returned.opts <- retval[["opts"]]
  predictions <- retval[["predictions"]]
  
  returned.opts$S_end <- 0
  returned.opts$E_end <- 0
  returned.opts$I_end <- 0
  returned.opts$R_end <- 0
  returned.opts$D_end <- 0
  for(i in 1:nrow(returned.opts)){
    p <- returned.opts$parameter.set[i]
    p <- predictions[[p]]
    p <- p[nrow(p),]
    returned.opts$S_end[i] <- p$S
    returned.opts$E_end[i] <- p$E
    returned.opts$I_end[i] <- p$I
    returned.opts$R_end[i] <- p$R
    returned.opts$D_end[i] <- p$D
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

slim.options <- function(final.opts){
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

nextstepSEIR <- function(simulation.start,simulation.end,last.step,R0,data){ 
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
    new.ops <- trySEIR(simulation.start,simulation.end,opts,data)
    for(j in 1:length(new.ops)){
      save.opts[[length(save.opts)+1]] <- last.step[[i]]
      save.opts[[length(save.opts)]] <- rbind(save.opts[[length(save.opts)]],new.ops[[j]])
    }
  }
  return(save.opts)
}

modelSEIR <- function(data,opts.R0,opts.E_0,opts.I_0){
names(data)[which(names(data)!="Date") & names(data)!="Underreporting")] <- "cases.observed"
data <- na.omit(data)
data$Date <- as.character(data$Date)
data$Date <- gsub("\\.","-",as.character(data$Date))
data$Date <- as.Date(data$Date,format="%d-%m-%Y")
data$day <- as.numeric(data$Date-min(data$Date))+1
data$days.to.today <- as.numeric(as.Date(format(Sys.time(), "%Y-%m-%d"))-data$Date)

today <- as.Date(format(Sys.time(), "%Y-%m-%d"))
days.to.today <- min(data$days.to.today)
days <- data$day[order(data$day)]
start.date <- max(days)-floor(max(days)/days.break)*days.break+1



opts <- expand.grid(opts.R0,opts.E_0,opts.I_0)
names(opts) <- c("R0","E_0","I_0")
opts$R_0 <- data$cases.observed[data$day==min(data$day[data$day>=start.date])]-opts$E_0
opts$D_0 <- 0

step1 <- trySEIR(start.date,start.date+days.break-1,opts,data)
step1 <- slim.options(step1)

first.day <- start.date+days.break
while(first.day<max(days)){
  print("--")
  print(first.day)
  print(max(days))
  print(length(step1))
  step1 <- nextstepSEIR(first.day,first.day+days.break-1,step1,opts.R0,data)
  step1 <- slim.options(step1)
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
                     D_0=old$D_end)
    sirmodel <- sirmodel[-1,]
    sirmodel$p <- old$p
    sirmodel$index <- old$index
    sirmodel$R0 <- old$R0
    
    results <- rbind(results,sirmodel)
  } 
    sirmodel <- SEIR(current$d.start-1,current$d.end,
                     R0=current$R0,
                     time.incubation=time.incubation,
                     time.infectedL=time.infectedL,
                     time.infectedD=time.infectedD,
                     cfr=cfr,
                     E_0=current$E_0,
                     I_0=current$I_0,
                     R_0=current$R_0,
                     D_0=current$D_0)
    sirmodel <- sirmodel[-1,]
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
                         cases.estimated2=sum(cases.estimated^2*adjp)),
                    by=list(time)]
collapsed$R0.sd <- sqrt(collapsed$R02-collapsed$R0^2)
collapsed$cases.estimated.sd <- sqrt(collapsed$cases.estimated2-collapsed$cases.estimated^2)
collapsed <- collapsed[,R02:=NULL]
collapsed <- collapsed[,cases.estimated2:=NULL]
setnames(collapsed,1,"day")

collapsed$cases.estimated.l95 <- collapsed$cases.estimated-1.96*collapsed$cases.estimated.sd
collapsed$cases.estimated.u95 <- collapsed$cases.estimated+1.96*collapsed$cases.estimated.sd

collapsed$R0.l95 <- collapsed$R0-1.96*collapsed$R0.sd
collapsed$R0.u95 <- collapsed$R0+1.96*collapsed$R0.sd

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

opts.R0 <- seq(0.2,4.0,0.2)
opts.E_0 <- seq(2,20,4)
opts.I_0 <- seq(2,20,4)

results.models <- list()
for(i in names(data.models)){
	results.models[[i]] <- modelSEIR(data.models[[i]],opts.R0,opts.E_0,opts.I_0)
}

try(results.models[["total"]]<-NULL,TRUE)
aggregated <- rbindlist(results.models)
aggregated <- aggregated[!is.na(aggregated$days.to.today),]
aggregated <- aggregated[,
	list(cases.estimated=sum(cases.estimated),
		cases.estimated.sd=sum(cases.estimated.sd^2),
		cases.observed=sum(cases.observed)),
	by=list(days.to.today)]



aggregated$cases.estimated.sd <- sqrt(aggregated$cases.estimated.sd)
aggregated$cases.estimated.l95 <- aggregated$cases.estimated-1.96*aggregated$cases.estimated.sd
aggregated$cases.estimated.u95 <- aggregated$cases.estimated+1.96*aggregated$cases.estimated.sd

results.models[["total"]] <- aggregated

#q <- ggplot(mapping=aes(x=day))
#q <- q + geom_line(data=collapsed,aes(y=R0),col="red")
#q

for(j in names(results.models)){
ybreaks <- c(50,100,1000,5000,10000,20000,40000,80000)
xbreaks <- seq(-28*3,28*2,28)
#text <- collapsed[nrow(collapsed),]
q <- ggplot(results.models[[j]],aes(x=days.to.today))
for(i in ybreaks) q <- q + geom_hline(yintercept=log2(i),lty=2,lwd=2)
for(i in xbreaks) q <- q + geom_vline(xintercept=i,lty=2,lwd=2)
q <- q + geom_ribbon(aes(ymin=log2(cases.estimated.l95),ymax=log2(cases.estimated.u95)),
                    fill="red",alpha=0.4)
q <- q + geom_point(aes(y=log2(cases.observed)),size=10)
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
ybreaks <- c(50,1000,5000,10000,20000,40000)
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





