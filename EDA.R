require(deSolve)
source("ggplot_theme_Publication.R")
options(scipen = 999)

DatasetCovid <-
  read.csv(
    url(
      'https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'
    )
  )

DatasetCovid <- DatasetCovid[0:650,]
Infetti <- DatasetCovid$totale_positivi
Rimossi <- DatasetCovid$dimessi_guariti + DatasetCovid$deceduti
date <- DatasetCovid$data

dati_reali <- data.frame(
  infetti = Infetti,
  rimossi = Rimossi
)

tempo <- as.Date(date)
colors <- c("Infetti" = 3)

ggplot(data = dati_reali, aes(x = tempo)) +
  geom_line(size = 0.8, aes(y = infetti, color = "Infetti")) +
  scale_color_manual("Dati",values = colors) +
  labs(title= "Confronto tra dati reali e modello SEIR stimato",
       subtitle=  "COVID-19 SEIR, Italia (2020/09/10 - 2021/09/05)",
       x="Tempo", y="Popolazione") +
  theme_Publication() +
  scale_x_date(breaks = date_breaks("months"))
