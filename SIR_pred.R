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

initDs <- 679 + 6
fineDs <- length(DatasetCovid[,1])
dsCovid <- DatasetCovid[initDs:fineDs,]

Rimossi <-
  dsCovid$dimessi_guariti + dsCovid$deceduti - dsCovid$dimessi_guariti[1] - dsCovid$deceduti[1]

Infetti <- dsCovid$totale_positivi

Pop <- 59000000
tempo <- 0:(length(Infetti) - 1)
NInit <- Pop
NrowLossArray <- 10
lossArray <- matrix(0, NrowLossArray, 4)

counter <- 1
exec_optim <- TRUE

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
  t <- tempo
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
  maxI <- max(out$I)
  diff <- sum(0.5 * ((out$I - Infetti)/maxI)^2 + 0.5 * ((out$R - Rimossi)/maxI)^2)
  sse <- diff
}

if (exec_optim) {
  init <- 0.1
  passo <- 0.01
  fine <- init + passo * (NrowLossArray - 1)
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
      lower = c(0, 0),
      upper = c(1, 1),
      hessian = TRUE
    )
    
    lossArray[counter, ] <- c(fit$par[1], fit$par[2], fit$value, prop)
    counter <- counter + 1
    cat("loading ", (counter - 1) * 100 /NrowLossArray,"%\n")
  }

  idxRes <- match(
    min(lossArray[,3]),
    lossArray[,3]
  )

  betaRes <- lossArray[idxRes,1]
  gammaRes <- lossArray[idxRes,2]
  N <- NInit * lossArray[idxRes,4]
  S0 <- N - Infetti[1] - Rimossi[1]

}else{
  # SIR Model
  N <- NInit * 0.109
  betaRes <- 0.1727
  gammaRes <- 0.0415582
  S0 <- N - Infetti[1] - Rimossi[1]
}

RZero <- betaRes/gammaRes

cat("il valore ottimo di R0 è: ", RZero,"\n")

Suscettibili <- N - Rimossi - Infetti

dati_reali <- data.frame(
  suscettibili = Suscettibili,
  infetti = Infetti,
  rimossi = Rimossi
)

I0 <- Infetti[1]
R0 <- Rimossi[1]

giorni_predizione <- 100
t <- seq(0, fineDs-initDs + giorni_predizione - 1, by = 1)

mod.pred <- as.data.frame(
  ode(
    c(S = S0, I = I0, R = R0),
    times = t,
    closed.sir.model,
    c(betaRes, gammaRes),
    hmax = 1 / 120
  )
)

res <- c(match( max(mod.pred$I), mod.pred$I) - length(Infetti), max(mod.pred$I))
cat("Il picco di infetti (", res[2],") si avrà in",res[1],"giorni\n")

colors <- c("Infetti" = 3,
            "Infetti SIR" = "green")

ggplot(data = dati_reali, mapping = aes(x = tempo, y = infetti)) +
  geom_point(aes(color = "Infetti")) +
  geom_line(data = mod.pred, mapping = aes(x = t, y = I, color = "Infetti SIR")) +
  scale_color_manual("Dati",values = colors) +
  labs(title= "Confronto tra Infetti reali e Infetti stimati",
       subtitle=  "COVID-19 Infetti, Italia (2020/10/10 - 2021/12/19)",
       x="Tempo", y="Infetti") +
  theme_Publication()

colors <- c("Suscettibili" = 2,
            "Infetti" = 3,
            "Rimossi" = 4,
            "Suscettibili SIR" = "red",
            "Infetti SIR" = "green",
            "Rimossi SIR" = "blue")

ggplot(data = dati_reali, aes(x = tempo)) +
  geom_point(size = 0.8, aes(y = suscettibili, color = "Suscettibili")) +
  geom_point(size = 0.8, aes(y = infetti, color = "Infetti")) +
  geom_point(size = 0.8, aes(y = rimossi, color = "Rimossi")) +
  geom_line(data = mod.pred, mapping = aes(x = t,y = S, color = "Suscettibili SIR")) +
  geom_line(data = mod.pred, mapping = aes(x = t,y = I, color = "Infetti SIR")) +
  geom_line(data = mod.pred, mapping = aes(x = t,y = R, color = "Rimossi SIR")) +
  scale_color_manual("Dati",values = colors) +
  labs(title= "Confronto tra dati reali e modello SIR stimato",
       subtitle=  "COVID-19 SIR, Italia (2020/10/10 - 2021/12/19)",
       x="Tempo", y="Popolazione") +
  theme_Publication()
