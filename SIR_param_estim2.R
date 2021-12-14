require(deSolve)

DatasetCovid <-
  read.csv(
    url(
      'https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'
    )
  )

Infected <-  DatasetCovid$totale_positivi[150:360]
Rimossi <-
  DatasetCovid$dimessi_guariti[150:360] + DatasetCovid$deceduti[150:360] - DatasetCovid$dimessi_guariti[150] - DatasetCovid$deceduti[150]

Day <- 0:(length(Infected) - 1)
N <- 59000000 * 0.05
lossArray <- matrix(0, 20, 4)

counter <- 1
exec_optim = FALSE

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
  cases <- Infected
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
  
  diff <- sum(abs((out$I - cases)))
  print(diff)
  sse <- diff
}

if (exec_optim) {
  for (prop in seq(1, 1.8, by = 0.1)) {
    N <- N * prop
    
    S0 <- N - Infected[1] - Rimossi[1]
    I0 <- Infected[1]
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
}

# SIR Model
N <- 59000000 * 0.05
beta <- 0.067
gamma <- 0.027

S0 <- N - Infected[1] - Rimossi[1]
I0 <- Infected[1]
R0 <- Rimossi[1]


plot(
  Day,
  Infected,
  type = 'b',
  xlab = 'Day',
  ylab = 'I(t)',
  col = 'red'
)

t <- seq(1, 360, by = 1)

mod.pred <- as.data.frame(ode(
  c(S = S0, I = I0, R = R0),
  times = t,
  closed.sir.model,
  c(beta, gamma),
  hmax = 1 / 120
))


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
  Infected,
  type = "b",
  pch = 15,
  add = TRUE,
  col = 5
)


legend(
  x = "right",
  y = 0.92,
  c("Susceptibles", "Infecteds", "Recovereds"),
  pch = 1,
  col = 2:4
)
