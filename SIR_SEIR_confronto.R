require(deSolve)
source("ggplot_theme_Publication.R")
options(scipen = 999)

DatasetCovid <-
  read.csv(
    url(
      'https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'
    )
  )


initDs <- 500
fineDs <- 610
dsCovid <- DatasetCovid[initDs:fineDs,]

Infetti <- dsCovid$totale_positivi
Rimossi <- dsCovid$dimessi_guariti + dsCovid$deceduti -
  dsCovid$dimessi_guariti[1] -
  dsCovid$deceduti[1]

Pop <- 59000000
tempo <- 0:(length(Infetti) - 1)
NInit <- Pop
NrowLossArray <- 10
lossArray <- matrix(0, NrowLossArray, 4)

counter <- 1

closed.sir.model <- function(t, x, params) {
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

closed.seir.model <- function (t, x, params) {
  S <- x[1]
  E <- x[2]
  I <- x[3]
  R <- x[4]

  beta <- params[1]
  gamma <- params[2]
  delta <- 1/7

  dS <- - (beta / N * S * I)
  dE <- (beta / N * S * I) - (delta * E)
  dI <- (delta * E) - (gamma * I)
  dR <- gamma * I
  list(c(dS, dE, dI, dR))
}

# SIR Model
N <- NInit * 0.1
betaRes <- 0.0844250
gammaRes <- 0.0373120
S0 <- N - Infetti[1] - Rimossi[1]

# SEIR Model
N_SEIR <- NInit * 0.1
betaRes_SEIR <- 0.115948
gammaRes_SEIR <- 0.0425646

Suscettibili <- N - Rimossi - Infetti

dati_reali <- data.frame(
  suscettibili = Suscettibili,
  infetti = Infetti,
  rimossi = Rimossi
)

I0 <- Infetti[1]
R0 <- Rimossi[1]

E0_SEIR <- Infetti[1] * 0.7
I0_SEIR <- Infetti[1] * 0.3
R0_SEIR <- Rimossi[1]


giorni_predizione <- 100
t <- seq(1, fineDs - initDs + giorni_predizione, by = 1)

mod.pred <- as.data.frame(
  ode(
    c(S = S0, I = I0, R = R0),
    times = t,
    closed.sir.model,
    c(betaRes, gammaRes),
    hmax = 1 / 120
  )
)

mod.pred_SEIR <- as.data.frame(
  ode(
    c(S = S0, E = E0_SEIR, I = I0_SEIR, R = R0_SEIR),
    times = t,
    closed.seir.model,
    c(betaRes_SEIR, gammaRes_SEIR),
    hmax = 1 / 120
  )
)

colors <- c("Suscettibili SIR" = "red",
            "Infetti SIR" = "limegreen",
            "Rimossi SIR" = "blue",
            "Suscettibili SEIR" = "purple",
            "Esposti SEIR" = "orange2",
            "Infetti SEIR" = "magenta",
            "Rimossi SEIR" = "cadetblue4"
)

ggplot(data = dati_reali, aes(x = tempo)) +
  geom_line(data = mod.pred, mapping = aes(x = t, y = S, color = "Suscettibili SIR")) +
  geom_line(data = mod.pred, mapping = aes(x = t, y = I, color = "Infetti SIR")) +
  geom_line(data = mod.pred, mapping = aes(x = t, y = R, color = "Rimossi SIR")) +
  geom_line(data = mod.pred_SEIR, mapping = aes(x = t,y = S, color = "Suscettibili SEIR")) +
  geom_line(data = mod.pred_SEIR, mapping = aes(x = t,y = E, color = "Esposti SEIR")) +
  geom_line(data = mod.pred_SEIR, mapping = aes(x = t,y = I, color = "Infetti SEIR")) +
  geom_line(data = mod.pred_SEIR, mapping = aes(x = t,y = R, color = "Rimossi SEIR")) +
  scale_color_manual("Dati", values = colors) +
  labs(title = "Confronto tra modello SIR e modello SEIR",
       subtitle = "COVID-19 SIR, Italia (2020/09/10 - 2021/09/05)",
       x = "Tempo", y = "Popolazione") +
  theme_Publication()
