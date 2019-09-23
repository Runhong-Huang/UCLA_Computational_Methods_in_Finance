library(ggplot2)



#Q2 

Q1 = c(101020, 101291,101024,100787,100611,100479,100380)

ggplot() + geom_point(aes(seq(0.3,0.9, length = 7),Q1,col = "")) +
  geom_line(aes(seq(0.3,0.9, length = 7),Q1,col = "")) +
  xlab("Keppa") + 
  ylab("MBS price")


Q2 = c(111989,109474,107189,105860,105570,100787,95578.7 )
ggplot() + geom_point(aes(seq(0.03,0.09, length = 7),Q2,col = "")) +
  geom_line(aes(seq(0.03,0.09, length = 7),Q2,col = "")) +
  xlab("r mean ") + 
  ylab("MBS price")

Q3 = c(100584, 100681,100787,100899,101018,101144,101276,101411,101549,101690,101831)
ggplot() + geom_point(aes(seq(0.1, 0.2 , by =  0.01),Q3,col = "")) +
  geom_line(aes(seq(0.1,0.2, by = 0.01),Q3,col = "")) +
  xlab("sigma ") + 
  ylab("MBS price")

