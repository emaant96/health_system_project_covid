require(deSolve)
require(ggplot2)
require(smooth)
source("ggplot_theme_Publication.R")
options(scipen = 999)

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

Infetti <- DatasetCovid$totale_positivi
VarInfetti <- DatasetCovid$variazione_totale_positivi
VarInfetti <- c(sma(VarInfetti,7)$fitted) * 5
Deceduti <- ave(DatasetCovid$deceduti, FUN=function(x) c(0,diff(x)))
Intensive <- DatasetCovid$terapia_intensiva
smaDec <- sma(Deceduti,7)
Deceduti <- c(smaDec$fitted)
Vaccinati <- DatasetVax$people_fully_vaccinated
Vaccinati <- c(integer(length(Infetti) - length(Vaccinati)), Vaccinati)

InfettiSc <- Infetti
DecedutiSc <- Deceduti * 1100
IntensiveSc <- Intensive * 210
VaccinatiSc <- Vaccinati*(max(Infetti)/59550000)
date <- DatasetCovid$data

dati_reali <- data.frame(
  infetti = Infetti,
  deceduti = DecedutiSc,
  intensive = IntensiveSc,
  vaccinati = VaccinatiSc,
  varinfetti = VarInfetti
)

tempo <- as.Date(date)
colors <- c("Infetti" = 3,
            "Deceduti" = 2,
            "Intensive" = 4,
            "Vaccinati"= 7,
            "Variazione Infetti" = 6)

ggplot(data = dati_reali, aes(x = tempo)) +
  geom_line(size = 0.8, aes(y = infetti, color = "Infetti")) +
  geom_line(size = 0.8, aes(y = deceduti, color = "Deceduti")) +
  geom_line(size = 0.8, aes(y = intensive, color = "Intensive")) +
  geom_line(size = 0.8, aes(y = vaccinati, color = "Vaccinati")) +
  geom_line(size = 0.8, aes(y = varinfetti, color = "Variazione Infetti")) +
  scale_color_manual("Dati",values = colors) +
  labs(title= "Andamento casi COVID-19",
       subtitle=  "Italia (2020 - 2021)",
       x="Data", y="Infetti") +
  theme_Publication() +
  scale_x_date(breaks = date_breaks("months"))+
  theme(axis.text.x = element_text(hjust= 1,angle=45))

cat("Infetti scala", max(InfettiSc)/max(Infetti), ":1\n")
cat("Deceduti scala", max(DecedutiSc)/max(Deceduti), ":1\n")
cat("Intensive scala", max(IntensiveSc)/max(Intensive), ":1\n")
cat("Vaccinati scala 1:", max(Vaccinati)/max(VaccinatiSc), "\n")