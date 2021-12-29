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

ggplot(data = dati_reali, aes(x = DatasetCovid$data, group = 1)) +
  geom_line(mapping = aes(y = infetti, color = "green")) +
  geom_line(mapping = aes(y = rimossi, color = "red")) +
  scale_x_discrete(breaks = unique(DatasetCovid$data)[seq(1,length(DatasetCovid$data),120)])
  labs(title = "COVID-19, Italia",
       # subtitle = "",
       x = "Tempo", y = "Popolazione") + theme_Publication() #  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))


