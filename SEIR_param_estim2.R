require(deSolve)
require(ggplot2)
source("ggplot_theme_Publication.R")
options(scipen=999)

DatasetCovid <-
  read.csv('./dpc-covid19-ita-andamento-nazionale.csv')


initDs <- 510
fineDs <- 605
dsCovid <- DatasetCovid[initDs:fineDs,]

Rimossi <-
  dsCovid$dimessi_guariti + dsCovid$deceduti - dsCovid$dimessi_guariti[1] - dsCovid$deceduti[1]

Infetti <- dsCovid$totale_positivi


Pop <- 59000000
tempo <- 0:(length(Infetti) - 1)
NInit <- Pop
NrowLossArray <- 10
lossArray <- matrix(0, NrowLossArray, 4)
delta <- 1/14

counter <- 1
exec_optim <- FALSE
#exec_optim <- TRUE

closed.seir.model <- function (t, x, params) {
  S <- x[1]
  E <- x[2]
  I <- x[3]
  R <- x[4]
  
  beta <- params[1]
  gamma <- params[2]
  delta <- delta
  #
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
  S0 <- S0
  E0 <- E0
  I0 <- I0
  R0 <- R0
  out <- as.data.frame(ode(
    y = c(S = S0, E = E0, I = I0, R = R0),
    times = t,
    closed.seir.model,
    parms = c(beta, gamma),
    hmax = 1 / 120
  ))
  
  diff <- sum(0.5 * (out$I - Infetti)^2 + 0.5 * (out$R - Rimossi)^2)
  print(diff)
  sse <- diff
}

if (exec_optim) {
  init <- 0.005
  passo <- 0.005
  fine <- init + passo * (NrowLossArray - 1)
  for (prop in seq(init, fine, by = passo)) {
    N <- NInit * prop
    
    S0 <- N - Infetti[1] - Infetti[1/delta] - Rimossi[1]
    E0 <- Infetti[1/delta]
    I0 <- Infetti[1]
    R0 <- Rimossi[1]
    
    params0 <- c(0.076, 0.027)
    
    fit <- optim(
      params0,
      sse.seir,
      method = "L-BFGS-B",
      lower = c(0, 0),
      upper = c(1, 1),
      hessian = TRUE
    )

    lossArray[counter, ] <- c(fit$par[1], fit$par[2], fit$value, prop)
    counter <- counter + 1
    cat("counter: ", counter)
  }

  idxRes <- match(
    min(lossArray[,3]),
    lossArray[,3]
  )

  betaRes <- lossArray[idxRes,1]
  gammaRes <- lossArray[idxRes,2]
  N <- NInit * lossArray[idxRes,4]
  S0 <- N - Infetti[1] - Infetti[1/delta] - Rimossi[1]

}else{
  # SEIR Model
  N <- NInit * 0.03
  betaRes <- 0.228
  gammaRes <- 0.026
  S0 <- N - Infetti[1] - Infetti[1/delta] - Rimossi[1]
}

RZero <- betaRes/gammaRes

cat("il valore ottimo di R0 è: ", RZero)

Suscettibili <- N - Rimossi - Infetti

dati_reali <- data.frame(
  suscettibili = Suscettibili,
  infetti = Infetti,
  rimossi = Rimossi
)

E0 <- Infetti[1/delta]
I0 <- Infetti[1]
R0 <- Rimossi[1]

giorni_predizione <- 100
t <- seq(0, fineDs-initDs + giorni_predizione-1, by = 1)

mod.pred <- as.data.frame(
  ode(
    c(S = S0, E = E0, I = I0, R = R0),
    times = t,
    closed.seir.model,
    c(betaRes, gammaRes),
    hmax = 1 / 120
  )
)

colors <- c("Infetti" = 3,
            "Infetti SEIR" = "green")

ggplot(data = dati_reali, mapping = aes(x = tempo, y = infetti)) +
  geom_point(aes(color = "Infetti")) +
  geom_line(data = mod.pred, mapping = aes(x = t, y = I, color = "Infetti SEIR")) +
  scale_color_manual("Dati",values = colors) +
  labs(title= "Confronto tra Infetti reali e Infetti stimati",
     subtitle=  "COVID-19 Infetti, Italia (2020/10/10 - 2020/12/19)",
     x="Tempo", y="Infetti") +
  theme_Publication()


colors <- c("Suscettibili" = 2,
            "Infetti" = 3,
            "Rimossi" = 4,
            "Suscettibili SEIR" = "red",
            "Esposti SEIR" = "orange2",
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
  labs(title= "Confronto tra dati reali e modello SEIR stimato",
       subtitle=  "COVID-19 SEIR, Italia (2020/10/10 - 2020/12/19)",
       x="Tempo", y="Popolazione") +
  theme_Publication()