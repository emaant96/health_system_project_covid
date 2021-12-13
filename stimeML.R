require(deSolve)

Infected <- read.csv(url('https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'))
Infected <-  Infected$totale_positivi[150:360]

dim <- 100
Day <- 0:(length(Infected)-1)
N <- 59000000
sir.model <- function (t, x, params) {
  S <- x[1]
  I <- x[2]
  R <- x[3]

  beta <- params[1]
  gamma <- params[2]

  dS <- -beta*S/N*I
  dI <- beta*S/N*I-gamma*I
  dR <- gamma*I
  list(c(dS,dI,dR))
}

t <- seq(1,max(Day),by=0.05)

results <- array(0,dim)
paramResults <- array(0,dim)

for( i in 1:dim){

  beta <- 10.227
  gamma <- 0.0053
  pars <- c(beta,gamma)

  if(FALSE){
    plot(Day,
       Infected,
       type='b',
       xlab='Day',
       ylab='I(t)',
       col='red'
    )
  }

  mod.pred<-as.data.frame(
    lsoda(
      c(S=max(Infected),I=1,R=0),
      times=Day,
      sir.model,
      pars)
  )
  lines(mod.pred$I~Day)
  paramResults[i] <- pars[1]/pars[2]
  results[i] <- sum(abs((mod.pred$I-Infected)))
}

plot(
  paramResults,
  results,
  pch = 19,
  col = 'red'
)

matplot(
  Day,
  mod.pred[,2:4],
  type = "l",
  xlab = "Time",
  ylab = "Susceptibles and Recovereds",
  main = "SIR Model, COVID-19, Italy, Estimated",
  lwd = 1,
  lty = 1,
  bty = "l",
  col = 2:4
)
legend(
  x= "right",
  y=0.92,
  c("Susceptibles", "Infecteds", "Recovereds"),
  pch = 1,
  col = 2:4
)