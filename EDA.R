require(deSolve)
source("ggplot_theme_Publication.R")
options(scipen = 999)

DatasetCovid <-
  read.csv('./dpc-covid19-ita-andamento-nazionale.csv')

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
  labs(title= "Andamento casi COVID-19",
       subtitle=  "Italia (2020 - 2021)",
       x="Data", y="Infetti") +
  theme_Publication() +
  scale_x_date(breaks = date_breaks("months"))+
  theme(axis.text.x = element_text(hjust= 1,angle=45))
