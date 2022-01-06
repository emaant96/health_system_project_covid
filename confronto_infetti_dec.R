require(deSolve)
require(ggplot2)
require(smooth)
source("test_miei/ggplot_theme_Publication.R")
options(scipen = 999)

DatasetCovid <-
  read.csv(
    url(
      'https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'
    )
  )

Infetti <- DatasetCovid$totale_positivi

Deceduti <- ave(DatasetCovid$deceduti, FUN=function(x) c(0,diff(x)))
Intensive <- DatasetCovid$terapia_intensiva
Infetti <- Infetti
Deceduti <- Deceduti * 1100
smaDec <- sma(Deceduti,7)
Deceduti <- c(smaDec$fitted)
Intensive <- Intensive * 230
date <- DatasetCovid$data

dati_reali <- data.frame(
  infetti = Infetti,
  deceduti = Deceduti,
  intensive = Intensive
)

tempo <- as.Date(date)
colors <- c("Infetti" = 3,"Deceduti" = 2, "Intensive" = 4)

ggplot(data = dati_reali, aes(x = tempo)) +
  geom_line(size = 0.8, aes(y = infetti, color = "Infetti")) +
  geom_line(size = 0.8, aes(y = deceduti, color = "Deceduti")) +
  geom_line(size = 0.8, aes(y = intensive, color = "Intensive")) +
  scale_color_manual("Dati",values = colors) +
  labs(title= "Andamento casi COVID-19",
       subtitle=  "Italia (2020 - 2021)",
       x="Data", y="Infetti") +
  theme_Publication() +
  scale_x_date(breaks = date_breaks("months"))+
  theme(axis.text.x = element_text(hjust= 1,angle=45))
