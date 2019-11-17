# 191117 - Stocks plot

setwd("~/Documents/GitHub/Random_scripts/Stocks")

install.packages("rjson")
library("rjson")
# https://financialmodelingprep.com/developer/docs/#Company-Financial-Ratios


symbol <- "TSLA"
symbol <- "AAPL"



chart <- fromJSON(file=paste("https://financialmodelingprep.com/api/v3/historical-price-full/", symbol, "?serietype=line", sep=""))

chart <- tail(unlist(chart[[2]]), 365*10)


chart <- data.frame("date"=as.Date(chart[seq(1, length(chart), 2)]), "price"=as.numeric(chart[seq(2, length(chart), 2)]))


plot(chart$date, chart$price, type='l', yaxt="n", ylab="US $", xlab="")
axis(2, las=1)

data <- fromJSON(file="https://financialmodelingprep.com/api/v3/financials/income-statement/TSLA?period=quarter")


financials <- unlist(data[[2]][1])

for (i in 2:length(data[[2]])){
financials <- rbind(financials, unlist(data[[2]][i]))
}
financials <- as.data.frame(financials)
financials$date <- as.Date(financials$date)

for(i in 2:ncol(financials)){
financials[,i] <- as.numeric(financials[,i])
}




names(financials)



barplot(rev(financials$Revenue/1000000))






