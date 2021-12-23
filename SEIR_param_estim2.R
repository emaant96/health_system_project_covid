require(deSolve)

DatasetCovid <-
  read.csv(
    url(
      'https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'
    )
  )


initDs <- 200
fineDs <- 360
dsCovid <- DatasetCovid[initDs:fineDs,]

Infetti <- dsCovid$totale_positivi
Rimossi <-
  dsCovid$dimessi_guariti + dsCovid$deceduti - dsCovid$dimessi_guariti[1] - dsCovid$deceduti[1]

Pop <- 59000000
tempo <- 0:(length(Infetti) - 1)
NInit <- Pop * 0.05
NrowLossArray <- 10
lossArray <- matrix(0, NrowLossArray, 5)

counter <- 1
#exec_optim <- FALSE
exec_optim <- TRUE

closed.seir.model <- function (t, x, params) {
  S <- x[1]
  E <- x[2]
  I <- x[3]
  R <- x[4]
  
  beta <- params[1]
  gamma <- params[2]
  delta <- params[3]
  
  dS <- - (beta / N * S * I)
  dE <- (beta / N * S * I) - (delta * E)
  dI <- (delta * E) - (gamma * I)
  dR <- gamma * I
  list(c(dS, dE, dI, dR))
}

sse.seir <- function(params0) {
  t <- tempo
  beta <- params0[1]
  gamma <- params0[2]
  delta <- params0[3]
  S0 <- S0
  E0 <- E0
  I0 <- I0
  R0 <- R0
  out <- as.data.frame(ode(
    y = c(S = S0, E = E0, I = I0, R = R0),
    times = t,
    closed.seir.model,
    parms = c(beta, gamma, delta),
    hmax = 1 / 120
  ))
  
  diff <- sum(0.5 * (out$I - Infetti)^2 + 0.25 * (out$R - Rimossi)^2 + 0.25 * (out$I + out$R - (Infetti + Rimossi))^2)
  print(diff)
  sse <- diff
}

if (exec_optim) {
  init <- 0.1
  passo <- 0.1
  fine <- init + passo * (NrowLossArray - 1)
  for (prop in seq(init, fine, by = passo)) {
    N <- NInit * prop
    
    S0 <- N - Infetti[1] - Rimossi[1]
    E0 <- Infetti[1] * 0.7
    I0 <- Infetti[1] * 0.3
    R0 <- Rimossi[1]
    
    params0 <- c(0.076, 0.027, 1/7)
    
    fit <- optim(
      params0,
      sse.seir,
      method = "L-BFGS-B",
      lower = c(1/1000, 1 / 42, 1/30),
      upper = c(1/5, 1 / 11, 1/5),
      hessian = TRUE
    )
    
    lossArray[counter, ] <- c(fit$par[1], fit$par[2], fit$par[3], fit$value, prop)
    counter <- counter + 1
    cat("counter: ", counter)
  }

  idxRes <- match(
    min(lossArray[,4]),
    lossArray[,4]
  )

  betaRes <- lossArray[idxRes,1]
  gammaRes <- lossArray[idxRes,2]
  deltaRes <- lossArray[idxRes,3]
  N <- NInit * lossArray[idxRes,5]
  S0 <- N - Infetti[1] - Rimossi[1]

}else{
  # SIR Model
  N <- Pop * 0.05 * 0.7
  betaRes <- 0.112
  gammaRes <- 0.03
  deltaRes <- 0.179
  S0 <- N - Infetti[1] - Rimossi[1]
}

RZero <- betaRes/gammaRes

cat("il valore ottimo di R0 è: ", RZero)

Suscettibili <- N - Rimossi - Infetti

dati_reali <- data.frame(
  suscettibili = Suscettibili,
  infetti = Infetti,
  rimossi = Rimossi
)

E0 <- Infetti[1] * 0.7
I0 <- Infetti[1] * 0.3
R0 <- Rimossi[1]
t <- seq(1, fineDs, by = 1)

mod.pred <- as.data.frame(
  ode(
    c(S = S0, E = E0, I = I0, R = R0),
    times = t,
    closed.seir.model,
    c(betaRes, gammaRes, deltaRes),
    hmax = 1 / 120
  )
)

ggplot(data = dati_reali, mapping = aes(x = tempo, y = infetti)) +
  geom_point(color = "red") +
  geom_line(data = mod.pred, mapping = aes(x = t, y = I))

colors <- c("Suscettibili" = 2,
            "Infetti" = 3,
            "Rimossi" = 4,
            "Suscettibili SEIR" = "red",
            "Esposti SEIR" = "yellow",
            "Infetti SEIR" = "green",
            "Rimossi SEIR" = "blue")

ggplot(data = dati_reali, aes(x = tempo)) +
  geom_point(size = 0.8, aes(y = suscettibili, color = "Suscettibili")) +
  geom_point(size = 0.8, aes(y = infetti, color = "Infetti")) +
  geom_point(size = 0.8, aes(y = rimossi, color = "Rimossi")) +
  geom_line(data = mod.pred, mapping = aes(x = t,y = S, color = "Suscettibili SEIR")) +
  geom_line(data = mod.pred, mapping = aes(x = t,y = E, color = "Esposti SEIR")) +
  geom_line(data = mod.pred, mapping = aes(x = t,y = I, color = "Infetti SEIR")) +
  geom_line(data = mod.pred, mapping = aes(x = t,y = R, color = "Rimossi SEIR")) +
  scale_color_manual("Dati",values = colors) +
  labs(title= "Confronto tra dati reali e modello SEIR stimato", x="Tempo", y="Popolazione")






if(FALSE){
matplot(
  t,
  mod.pred[, 2:4],
  type = "l",
  xlab = "Tempo",
  ylab = "Popolazione",
  main = "COVID-19 SEIR, Italia (2020-09-10 - 2021-09-05)",
  lwd = 1,
  lty = 1,
  bty = "l",
  col = 2:4
)

matplot(
  tempo,
  dati_reali,
  type = "b",
  pch = 15,
  add = TRUE,
  col = 2:4
)

legend(
  x = "right",
  y = 0.92,
  c("Suscettibili", "Infetti", "Rimossi"),
  pch = 1,
  col = 2:4
)
}