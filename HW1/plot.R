library(readr)
library(ggplot2)


#Q2 
Q2 = read_csv("Q2.csv")

hist(as.numeric(unlist(Q2)), xlab = "X", main = "Histgram of X", col = 'blue')

ggplot() + geom_bar(aes(as.numeric(unlist(Q2)))) + xlab("X") + title("Q2 Distribution")

#Q3 
Q3 = read_csv("Q3.csv")

ggplot() + geom_bar(aes(as.numeric(unlist(Q3)))) + xlab("X")


#Q4
Q4 = read_csv("Q4.csv")

ggplot() + geom_histogram(aes(as.numeric(unlist(Q4))), bins = 100) + xlab("X")



       