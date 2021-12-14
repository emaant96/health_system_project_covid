require(deSolve)

DatasetCovid <-
  read.csv(
    url(
      'https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'
    )
  )

Infetti <- DatasetCovid$totale_positivi[150:360]
Rimossi <-
  DatasetCovid$dimessi_guariti[150:360] + DatasetCovid$deceduti[150:360] - DatasetCovid$dimessi_guariti[150] - DatasetCovid$deceduti[150]

Pop <- 59000000
Day <- 0:(length(Infetti) - 1)
NInit <- Pop * 0.05
NrowLossArray <- 8
lossArray <- matrix(0, NrowLossArray, 4)

counter <- 1
exec_optim <- FALSE

closed.sir.model <- function (t, x, params) {
  S <- x[1]
  I <- x[2]
  R <- x[3]
  
  beta <- params[1]
  gamma <- params[2]
  
  dS <- -beta / N * S * I
  dI <- beta / N * S * I - gamma * I
  dR <- gamma * I
  list(c(dS, dI, dR))
}

sse.sir <- function(params0) {
  t <- Day
  cases <- Infetti
  beta <- params0[1]
  gamma <- params0[2]
  S0 <- S0
  I0 <- I0
  R0 <- R0
  out <- as.data.frame(ode(
    y = c(S = S0, I = I0, R = R0),
    times = t,
    closed.sir.model,
    parms = c(beta, gamma),
    hmax = 1 / 120
  ))
  
  diff <- sum((out$I - cases)^2)
  print(diff)
  sse <- diff
}

if (exec_optim) {
  init <- 1
  passo <- 0.1
  fine <- init + passo * NrowLossArray;
  for (prop in seq(init, fine, by = passo)) {
    N <- NInit * prop
    
    S0 <- N - Infetti[1] - Rimossi[1]
    I0 <- Infetti[1]
    R0 <- Rimossi[1]
    
    params0 <- c(0.076, 0.027)
    
    fit <- optim(
      params0,
      sse.sir,
      method = "L-BFGS-B",
      lower = c(0.001, 1 / 42),
      upper = c(1, 1 / 11)
    )

    print(fit$par)
    
    lossArray[counter, ] <- c(fit$par[1], fit$par[2], fit$value, prop)
    counter <- counter + 1
    cat("counter: ", counter)
  }

  idxRes <- match(min(lossArray[,3]),lossArray[,3])

  betaRes <- lossArray[idxRes,1]
  gammaRes <- lossArray[idxRes,2]
  N <- NInit * lossArray[idxRes,4]
  S0 <- N - Infetti[1] - Rimossi[1]

}else{
  # SIR Model
  N <- Pop * 0.05 * 1.3
  betaRes <- 0.074# 0.067
  gammaRes <- 0.035# 0.027
  S0 <- N - Infetti[1] - Rimossi[1]
}

Suscettibili <- N - Rimossi - Infetti
dati_reali <- cbind(cbind(Suscettibili,Infetti),Rimossi)

I0 <- Infetti[1]
R0 <- Rimossi[1]

plot(
  Day,
  Infetti,
  type = 'b',
  xlab = 'Day',
  ylab = 'I(t)',
  col = 'red'
)

t <- seq(1, 360, by = 1)

mod.pred <- as.data.frame(
  ode(
    c(S = S0, I = I0, R = R0),
    times = t,
    closed.sir.model,
    c(betaRes, gammaRes),
    hmax = 1 / 120
  )
)

lines(mod.pred$I ~ t)

matplot(
  t,
  mod.pred[, 2:4],
  type = "l",
  xlab = "Time",
  ylab = "Susceptibles and Recovereds",
  main = "COVID-19 SIR Model, Italy (2020-07-22 - 2021-07-22)",
  lwd = 1,
  lty = 1,
  bty = "l",
  col = 2:4
)

matplot(
  Day,
  dati_reali,
  type = "b",
  pch = 15,
  add = TRUE,
  col = 2:4
)

legend(
  x = "right",
  y = 0.92,
  c("Susceptibles", "Infetti", "Recovereds"),
  pch = 1,
  col = 2:4
)
