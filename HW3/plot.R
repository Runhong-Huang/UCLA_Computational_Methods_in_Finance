library(ggplot2)

Q31 = c(0.08449,0.190743,0.380164,0.671287,1.07572,1.60126,2.23838,2.97144,3.77936,4.64575,5.5542)

ggplot() + geom_point(aes(x=seq(15:25), y = Q31)) + 
  xlab("S0") + 
  ylab("Call option Price")



Q32 = c(0.0857524,0.194263,0.382816,0.673219,1.07868,1.6016,2.23451,2.96303,3.76954,4.63619,5.54706)
ggplot() + geom_point(aes(x=seq(15:25), y = Q32)) + 
  xlab("S0") + 
  ylab("Call option Price")

D =read_csv("D.csv", col_names = F)
T =read_csv("T.csv", col_names = F)
V =read_csv("V.csv", col_names = F)
G =read_csv("G.csv", col_names = F)
R =read_csv("R.csv", col_names = F)

ggplot() + geom_point(aes(x=c(15:25), y = V$X1)) + 
  xlab("S0") + 
  ylab("Vega")


Q51 = read_csv("Q527.csv",col_names = F)
ggplot() + geom_point(aes(x=Q51$X1, y = Q51$X2)) + 
  xlab("X") + 
  ylab("Y")


