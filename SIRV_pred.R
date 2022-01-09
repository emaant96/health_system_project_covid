require(deSolve)
library(ggplot2)
source("ggplot_theme_Publication.R")
options(scipen = 999)

# 9.52, 9.08, 8.33, 8.78, 9.08

DatasetCovid <-
  read.csv(
    url(
      'https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'
    )
  )

DatasetVax <-
  read.csv(
    url(
      'https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/country_data/Italy.csv'
    )
  )

initDs <- 673
fineDs <- length(DatasetCovid[,1])
dsCovid <- DatasetCovid[initDs:fineDs,]

Rimossi <-
  dsCovid$dimessi_guariti + dsCovid$deceduti - dsCovid$dimessi_guariti[1] - dsCovid$deceduti[1]

Infetti <- dsCovid$totale_positivi
Vaccinati <- DatasetVax$total_boosters
Vaccinati <- Vaccinati[366:length(Vaccinati)]

Pop <- 59550000
tempo <- 0:(length(Infetti) - 1)
NInit <- Pop
NrowLossArray <- 10
lossArray <- matrix(0, NrowLossArray, 5)

counter <- 1
exec_optim <- TRUE

closed.sirv.model <- function (t, x, params) {
  S <- x[1]
  I <- x[2]
  R <- x[3]
  V <- x[4]
  
  beta <- params[1]
  gamma <- params[2]
  mu <- params[3]

  dS <- -beta / N * S * I - mu * S
  dI <- beta / N * S * I - gamma * I
  dR <- gamma * I
  dV <- mu * S
  list(c(dS, dI, dR, dV))
}

sse.sirv <- function(params0) {
  t <- tempo
  beta <- params0[1]
  gamma <- params0[2]
  mu <- params0[3]

  S0 <- S0
  I0 <- I0
  R0 <- R0
  V0 <- V0

  out <- as.data.frame(ode(
    y = c(S = S0, I = I0, R = R0, V = V0),
    times = t,
    closed.sirv.model,
    parms = c(beta, gamma, mu),
    hmax = 1 / 120
  ))

  diff <- sum(((out$I - Infetti)/S0)^2 + ((out$R - Rimossi)/S0)^2 + ((out$V - Vaccinati)/S0)^2)
  sse <- diff
}

if (exec_optim) {
  init <- 0.1
  passo <- 0.1
  fine <- init + passo * (NrowLossArray - 1)
  for (prop in seq(init, fine, by = passo)) {
    N <- NInit * prop
    
    S0 <- N - Infetti[1] - Rimossi[1] - Vaccinati[1]
    I0 <- Infetti[1]
    R0 <- Rimossi[1]
    V0 <- Vaccinati[1]
    
    params0 <- c(0.076, 0.027, 0.5)
    
    fit <- optim(
      params0,
      sse.sirv,
      method = "L-BFGS-B",
      lower = c(0, 0, 0),
      upper = c(1, 1, 1),
      hessian = TRUE
    )
    
    lossArray[counter, ] <- c(fit$par[1], fit$par[2], fit$par[3], fit$value, prop)
    counter <- counter + 1
    cat("loading ", (counter - 1) * 100 /NrowLossArray,"%\n")
  }

  idxRes <- match(
    min(lossArray[,4]),
    lossArray[,4]
  )

  betaRes <- lossArray[idxRes,1]
  gammaRes <- lossArray[idxRes,2]
  muRes <- lossArray[idxRes,3]
  N <- NInit * lossArray[idxRes,5]
  S0 <- N - Infetti[1] - Rimossi[1] - Vaccinati[1]

}else{
  # SIR Model
  N <- NInit * 0.06
  betaRes <- 0.21
  gammaRes <- 0.037
  S0 <- N - Infetti[1] - Rimossi[1]
}

RZero <- betaRes/gammaRes

cat("il valore ottimo di R0 è: ", RZero,"\n")

Suscettibili <- N - Rimossi - Infetti - Vaccinati

dati_reali <- data.frame(
  suscettibili = Suscettibili,
  infetti = Infetti,
  rimossi = Rimossi,
  vaccinati = Vaccinati
)

I0 <- Infetti[1]
R0 <- Rimossi[1]
V0 <- Vaccinati[1]

giorni_predizione <- 100
t <- seq(0, fineDs-initDs + giorni_predizione - 1, by = 1)

mod.pred <- as.data.frame(
  ode(
    c(S = S0, I = I0, R = R0,V = V0),
    times = t,
    closed.sirv.model,
    c(betaRes, gammaRes, muRes),
    hmax = 1 / 120
  )
)

res <- c(match( max(mod.pred$I), mod.pred$I) - length(Infetti), max(mod.pred$I))
cat("Il picco di infetti (", res[2],") si avrà in",res[1],"giorni\n")

colors <- c("Suscettibili" = 2,
            "Infetti" = 3,
            "Rimossi" = 4,
            "Vaccinati" = 7,
            "Suscettibili SIRV" = "red",
            "Infetti SIRV" = "green",
            "Rimossi SIRV" = "blue",
            "Vaccinati SIRV" = "yellow")

ggplot(data = dati_reali, aes(x = tempo)) +
  geom_point(size = 0.8, aes(y = suscettibili, color = "Suscettibili")) +
  geom_point(size = 0.8, aes(y = infetti, color = "Infetti")) +
  geom_point(size = 0.8, aes(y = rimossi, color = "Rimossi")) +
  geom_point(size = 0.8, aes(y = vaccinati, color = "Vaccinati")) +
  geom_line(data = mod.pred, mapping = aes(x = t,y = S, color = "Suscettibili SIRV")) +
  geom_line(data = mod.pred, mapping = aes(x = t,y = I, color = "Infetti SIRV")) +
  geom_line(data = mod.pred, mapping = aes(x = t,y = R, color = "Rimossi SIRV")) +
  geom_line(data = mod.pred, mapping = aes(x = t,y = V, color = "Vaccinati SIRV")) +
  scale_color_manual("Dati",values = colors) +
  labs(title= "Confronto tra dati reali e modello SIRV stimato",
       subtitle=  "COVID-19 SIRV, Italia (2020/10/10 - 2021/12/19)",
       x="Tempo", y="Popolazione") +
  theme_Publication()
