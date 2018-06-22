library(ggplot2)
library(forecast)
library(tseries)
library(fBasics)
library(lmtest)
library(zoo)
library(fUnitRoots)

#Read data
ds <- read.table("D:/CSC 425/seaice.csv", header = T, sep = ',')
df <- subset(ds, select = -c(6))

df$hemisphere <- gsub('north', 'N', df$hemisphere)
df$hemisphere <- gsub('south', 'S', df$hemisphere)

df$Date <- as.Date(with(df, paste(Year, Month, Day, sep = '-')), "%Y-%m-%d")

#Split into two dataframe based on hemisphere N or S
mydf <- subset(df, select = -c(1:3))
mydf <- mydf[, c(4, 1, 2, 3)]
mydfN <- split(mydf, mydf$hemisphere)[['N']]
mydfS <- split(mydf, mydf$hemisphere)[['S']]

newDate <- mydfN$Date
newExtent <- (mydfN$Extent + mydfS$Extent)/2

myd <- data.frame('Date' = newDate, 'Extent' = newExtent)

extts = ts(myd[, 2], start = c(1978, 10), end  = c(2017, 6), freq = 12)
extAts = ts(myd[, 2], start = c(1978, 10), end  = c(2017, 6), freq = 1)
basicStats(myd$Extent)

#Histogram and Q-Q plot
hist(myd$Extent, xlab = "Extent change", prob = TRUE, main = "Histogram")
xfit <- seq(min(myd$Extent), max(myd$Extent), length = 40)
yfit <- dnorm(xfit, mean = mean(myd$Extent), sd = sd(myd$Extent))
lines(xfit, yfit, col="blue", lwd = 2)

qqnorm(myd$Extent)
qqline(myd$Extent, col = 2)

plot(extts, ylab = 'Extent change')

acf(myd$Extent)
pacf(myd$Extent)

#First difference
adfTest(myd$Extent, lags = 5, type = c("ct"))
diffvar = diff(extts)
plot(diffvar)
adfTest(coredata(diffvar), lags = 5, type=c("c"))


hist(diffvar, xlab = "Extent change", prob = TRUE, main = "Histogram")
xfit <- seq(min(diffvar), max(diffvar), length = 40)
yfit <- dnorm(xfit, mean = mean(diffvar), sd = sd(diffvar))
lines(xfit, yfit, col="blue", lwd = 2)

qqnorm(diffvar)
qqline(diffvar, col = 2)

#ARMIA Model
auto.arima(extts, ic =c("bic"), trace = TRUE, allowdrift = TRUE)
m1 = Arima(extts, order = c(2,1,1), method = 'ML')
coeftest(m1)

acf(m1$residuals)

Box.test(m1$residuals,lag = 3,type = 'Ljung-Box', fitdf = 1)
Box.test(m1$residuals,lag = 6,type = 'Ljung-Box', fitdf = 1)

f = forecast(m1, h = 12)
f
plot(f, include = 200)

#SARIMA
x <- myd$Extent
ets = ts(x, frequency = 12, start = c(1978, 10), end  = c(2017, 6))
plot(ets, type = 'l')

hist(x, xlab = "Extent", freq = F)
xfit <- seq(min(x), max(x), length = 40)
yfit <- dnorm(xfit, mean = mean(x), sd = sd(x))
lines(xfit, yfit, col = "black", lwd = 2)

par(mfcol = c(1, 1))
acf(as.vector(ets), lag.max = 30, main = "ACF")

dx = diff(ets)
acf(as.vector(dx), lag.max = 26, main = "ACF of DX starts")

sdx = diff(dx, 12)
acf(as.vector(sdx), lag.max = 30, main = "ACF of DSDX starts")
pacf(as.vector(sdx), lag.max = 30, main = "PACF of DSDX starts")

m2 = Arima(ets, order = c(2, 1, 1), seasonal = list(order = c(1, 1, 0), period = 12), method = "ML")
m2
coeftest(m1)
acf(m1$residuals)

Box.test(m1$residuals, 4, "Ljung-Box", fitdf = 2)
Box.test(m1$residuals, 12, "Ljung-Box", fitdf = 2) 

f1 = forecast(m2, h = 12)
f1
plot(f1, include = 200)
lines(ts(c(f1$fitted, f1$mean), frequency = 12, start = c(1978, 10), end = c(2017, 6)), col = "blue")

#Backtest
source('D:/CSC 425/backtest.R')
pm1 = backtest(m1, extts, 200, 1)
pm2 = backtest(m2, extts, 200, 1)

#Annual Change (ARIMA)
auto.arima(extAts, ic =c("aic"), trace = TRUE, allowdrift = TRUE)
mA1 = Arima(extAts, order = c(0,2,1), method = 'ML')
coeftest(mA1)

Box.test(m1$residuals,lag = 3,type = 'Ljung-Box', fitdf = 1)
Box.test(m1$residuals,lag = 6,type = 'Ljung-Box', fitdf = 1)

fA = forecast(mA1, h = 5)
plot(fA, include = 200)

#Tentative SARIMA Model (North pole)
xN <- mydfN$Extent
nts = ts(xN, frequency = 12, start = c(1978, 10), end  = c(2017, 6))
plot(nts, type = 'l')

hist(xN, xlab = "Extent", freq = F)
xfit <- seq(min(xN), max(xN), length = 40)
yfit <- dnorm(xfit, mean = mean(xN), sd = sd(xN))
lines(xfit, yfit, col = "black", lwd = 2)

par(mfcol = c(1, 1))
acf(as.vector(nts), lag.max = 30, main = "ACF")

ndx = diff(nts)
acf(as.vector(ndx), lag.max = 26, main = "ACF of DX extents")

nsdx = diff(ndx, 12)
acf(as.vector(nsdx), lag.max = 30, main = "ACF of DSDX extents")

mn = Arima(nts, order = c(2, 0, 0), seasonal = list(order = c(0, 1, 2), period = 4), method = "ML")
coeftest(mn)
acf(mn$residuals)

Box.test(mn$residuals, 4, "Ljung-Box", fitdf = 2)
Box.test(mn$residuals, 12, "Ljung-Box", fitdf = 2) 

fn = forecast(mn, h = 10)
plot(fn, include = 1000)
lines(ts(c(fn$fitted, fn$mean), frequency = 12, start = c(1978, 10), end = c(2017, 6)), col = "blue")

#Tentative SARIMA Model (South pole)
xS <- mydfS$Extent
sts = ts(xS, frequency = 12, start = c(1978, 10), end = c(2017, 6))
plot(sts, type = 'l')

par(mfcol = c(1, 1))
acf(as.vector(sts), lag.max = 30, main = "ACF")

sdx = diff(sts)
acf(as.vector(sdx), lag.max = 26, main = "ACF of DX extents")

ssdx = diff(sdx, 12)
acf(as.vector(ssdx), lag.max = 30, main = "ACF of DSDX extents")

ms = Arima(sts, order = c(1, 0, 0), seasonal = list(order = c(2, 1, 0), period = 4), method = "ML")
coeftest(ms)
acf(ms$resid)

Box.test(ms$residuals, 4, "Ljung-Box", fitdf = 3)
Box.test(ms$residuals, 12, "Ljung-Box", fitdf = 3) 

fs = forecast(ms, h = 10)
plot(fs, include = 1000)
lines(ts(c(fs$fitted, fs$mean), frequency = 12, start = c(1978, 10), end = c(2017, 6)), col = "blue")
