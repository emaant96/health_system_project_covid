require(deSolve)
source("ggplot_theme_Publication.R")
options(scipen = 999)

DatasetCovid <-
  read.csv(
    url(
      'https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'
    )
  )

Infetti <- DatasetCovid$totale_positivi
Rimossi <- DatasetCovid$dimessi_guariti + DatasetCovid$deceduti

dati_reali <- data.frame(
  infetti = Infetti,
  rimossi = Rimossi
)

ggplot(data = dati_reali, aes(x = tempo)) +
  geom_point(size = 0.8, aes(y = infetti, color = "Infetti"))+
  scale_color_manual("Dati",values = colors) +
  labs(title= "Confronto tra dati reali e modello SEIR stimato",
       subtitle=  "COVID-19 SEIR, Italia (2020/09/10 - 2021/09/05)",
       x="Tempo", y="Popolazione") +
  theme_Publication()

