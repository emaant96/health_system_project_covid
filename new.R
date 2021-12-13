####
####
####

library(deSolve)

Infected <- read.csv(url('https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'))
Infected <- Infected$totale_positivi[1:150]
day <- 0:(length(Infected)-1)
N <- 59000000

init <- c(S = max(Infected), I = Infected[1], R = 0)
plot(day, Infected)

SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta * (S/N) * I
    dI <- beta * (S/N) * I - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

RSS.SIR <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = day, func = SIR, parms = parameters)
  fit <- out[ , 3]
  RSS <- sum((Infected - fit)^2)
  return(RSS)
}

fit <- optim(
  c(0.9, 0.9), RSS.SIR, method = "BFGS"
)

fit$par

Opt <- setNames(fit$par, c("beta", "gamma"))

parameters <- c(beta = fit$par[1], gamma = fit$par[2])
parameters
times <- seq(0, length(Infected), by = 1)
out <- as.data.frame(ode(y = init, times = times,
                         func = SIR, parms = parameters))
out$time <- NULL

## PLOT SIR Data and Model
par(mfrow=c(2,1))

matplot(day,
        format(Infected / N, scientific = F, digits = 3),
        type = "l",
        xlab = "Time",
        ylab = "Infected",
        main = "Infected, COVID-19, Italy",
        lwd = 1,
        lty = 1,
        bty = "l",
        col = 3
)

legend(x= "bottomright",
       y=0.92, "Infecteds",
       pch = 1,
       col = 3
)

matplot(times,
        out,
        type = "l",
        xlab = "Time",
        ylab = "Susceptibles and Recovereds",
        main = "SIR Model, COVID-19, Italy, Estimated",
        lwd = 1,
        lty = 1,
        bty = "l",
        col = 2:4
)

legend(x= "right",
       y=0.92,
       c("Susceptibles", "Infecteds", "Recovereds"),
       pch = 1,
       col = 2:4
)