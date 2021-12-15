require(deSolve)
if (!require("forecast"))
  install.packages("forecast")
library(forecast)

DatasetCovid <-
  read.csv(
    url(
      'https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'
    )
  )

initDs <- 200
fineDs <- 360
Infetti <- DatasetCovid$totale_positivi[initDs:fineDs]
Rimossi <-
  DatasetCovid$dimessi_guariti[initDs:fineDs] + DatasetCovid$deceduti[initDs:fineDs] - DatasetCovid$dimessi_guariti[initDs] - DatasetCovid$deceduti[initDs]


## partition into train and test
train_series = Infetti[1:140]
test_series = Infetti[141:160]

AutoArimaModel = auto.arima(train_series)
print(AutoArimaModel)

# accmeasures=regr.eval(test_series, forecast$pred)


# forecast next 200 days
forecast = forecast(AutoArimaModel, h = 20)
#plot(forecast, shadecol=rgb(0,0,1,.4), flwd=1,
#     main="Forecasts of Covid-19",
#     xlab="Days", ylab="# Infected")
test_series_padding <- rep(NA, length(train_series));
test_series <- c(test_series_padding, test_series);

print(autoplot(forecast) + autolayer(ts(test_series)))

if (FALSE) {
matplot(
  seq(1,80, by = 1),
  forecast$mean,
  type = "l",
  xlab = "Time",
  ylab = "Susceptibles and Recovereds",
  main = "COVID-19 SIR Model, Italy (2020-09-10 - 2021-09-05)",
  lwd = 1,
  lty = 1,
  bty = "l",
  col = 2:4
)

matplot(
  seq(1,160, by = 1),
  test_series,
  type = "b",
  pch = 15,
  # add = TRUE,
  col = 2:4
)

legend(
  x = "right",
  y = 0.92,
  c("Susceptibles", "Infetti", "Recovereds"),
  pch = 1,
  col = 2:4
)
}

