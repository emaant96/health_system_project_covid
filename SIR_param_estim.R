# SIR model parameter estimation in R
## we want to estimate beta and gamma and then solve the ode with these values.

Infected <- read.csv(url('https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'))
Infected <- Infected$totale_positivi[1:150]

Day <- 0:(length(Infected)-1)
N <- 59000000 #pop of china

SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

library(deSolve)
init <- c(S = max(Infected), I = Infected[1], R = 0)

## LOSS FUNCTION
RSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = Day, func = SIR, parms = parameters)
  fit <- out[ , 2]
  sum((Infected - fit)^2)
}

Opt <- optim(c(0.9, 0.9), RSS, method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1))

Opt_par <- setNames(Opt$par, c("beta", "gamma"))
Opt_par # beta = 0.6134688, gamma = 0.3865312

# init <- c(S = 1-1e-6, I = 1e-6, 0.0)
parameters <- c(beta = Opt[["par"]][1], gamma = Opt[["par"]][2])
times <- seq(0, length(Infected), by = 1)
out <- as.data.frame(ode(y = init, times = times,
                         func = SIR, parms = parameters))
out$time <- NULL

## PLOT SIR Data and Model
par(mfrow=c(2,1))

matplot(Day, format(Infected / N, scientific = F, digits = 3),
        type = "l", xlab = "Time", ylab = "Infected",
        main = "Infected, COVID-19, China", lwd = 1, lty = 1, bty = "l", col = 3)
legend(x= "bottomright", y=0.92, "Infecteds", pch = 1, col = 3)

matplot(times, out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds",
        main = "SIR Model, COVID-19, China, Estimated", lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(x= "right", y=0.92, c("Susceptibles", "Infecteds", "Recovereds"),
       pch = 1, col = 2:4)
