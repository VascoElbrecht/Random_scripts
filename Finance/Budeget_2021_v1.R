


data <- read.csv("~/Downloads/2020 portfolio - Expenses (3).csv")



head(data)

data <- data[data$Type=="Groceries",]



x <- aggregate(data$CHF, list(data$Week), "sum")

m <- barplot(rev(x[,2]*-1), ylab="CHF")
axis(1, m, labels=x[,1], tick=F)









