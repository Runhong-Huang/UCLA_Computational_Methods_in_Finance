library(ggplot2)

Q31 = c(0.08449,0.190743,0.380164,0.671287,1.07572,1.60126,2.23838,2.97144,3.77936,4.64575,5.5542)

ggplot() + geom_point(aes(x=seq(15:25), y = Q31)) + 
  xlab("S0") + 
  ylab("Call option Price")


