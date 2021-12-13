Infected <- read.csv(url('https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'));
Infected <- Infected$totale_positivi[225:350]

Day <- 0:(length(Infected)-1)
N <- 59000000 #pop of china

closed.sir.model <- function (t, x, params) {
  S <- x[1]
  I <- x[2]
  R <- x[3]

  beta <- params[1]
  gamma <- params[2]
  dS <- -beta/N*S*I
  dI <- beta/N*S*I-gamma*I
  dR <- gamma*I
  list(c(dS,dI,dR))
}

require(deSolve)
sse.sir <- function(params0){
  t <- Day
  cases <- Infected
  beta <- params0[1]
  gamma <- params0[2]
  S0 <- max(Infected)
  I0 <- 1
  R0 <- 0
  out <- as.data.frame(ode(y=c(S=S0,I=I0,R=R0),times=t,closed.sir.model,parms=c(beta,gamma),hmax=1/120))
  sse<-sum((out$I-cases)^2)
}

params0<-c(0.484252774, 0.005837638)
fit0 <- optim(params0,sse.sir, method = "BFGS"); fit0$par

plot(Day,Infected, type='b', xlab='Day', ylab='I(t)',col='red')
t <- seq(1,max(Day),by=0.05)
mod.pred<-as.data.frame(lsoda(c(S=max(Infected),I=1,R=0),times=t,closed.sir.model,fit0$par,hmax=1/120))
lines(mod.pred$I~t)

out <- as.data.frame(ode(y=c(S=max(Infected),I=1,R=0),times=t,closed.sir.model,parms=c(35.37, 0.0058)))

matplot(t, out[,2:4], type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "SIR Model, COVID-19, China, Estimated", lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(x= "right", y=0, c("Susceptibles", "Infecteds", "Recovereds"), pch = 1, col = 2:4)
