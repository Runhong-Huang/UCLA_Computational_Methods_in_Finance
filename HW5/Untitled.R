library(ggplot2)



#Q2 
Q1 = c (1877.52,2710.64,3524.73,4228.75,4833.31,5341.85,2823.76,3631.79,4296.55,4874.73,5229.68,5680.84,3416.42,4011.11,4661.04,5074.21,5379.15,5633.39,3827.14,4277.26,4753.47,5081.28,5305.04,5528.77,3970.62,4362.9,4733.53,4989.49,5211.04,5442.81,4071.09,4358.72,4679.96,4911.96,5144.16,5294.6,4093.7,4355.34,4624.25,4869.3,5094.92,5199.65,4013.62,4278.99,4527.28,4761.13,4975.64,5136.75,4013.25,4200.52,4480.13,4698.97,4919.03,5075.87)

T = c(seq(3,8 ,by = 1))

k = matrix(0, nrow = 9, ncol = 6)

for (i in seq(1,54,6)){
  print(i)
  k[(i-1)/6 +1,] = Q1[i:(i+5)]
}

ggplot() + geom_point(aes(T,k[1,],col = "0")) +geom_line(aes(T,k[1,],col = "0")) +
  geom_point(aes(T,k[2,],col = "0.1")) +geom_line(aes(T,k[2,],col = "0.1")) +
  geom_point(aes(T,k[3,],col = "0.2")) +geom_line(aes(T,k[3,],col = "0.2")) +
  geom_point(aes(T,k[4,],col = "0.3")) +geom_line(aes(T,k[4,],col = "0.3")) +
  geom_point(aes(T,k[5,],col = "0.4")) +geom_line(aes(T,k[5,],col = "0.4")) +
  geom_point(aes(T,k[6,],col = "0.5")) +geom_line(aes(T,k[6,],col = "0.5")) +
  geom_point(aes(T,k[7,],col = "0.6")) +geom_line(aes(T,k[7,],col = "0.6")) +
  geom_point(aes(T,k[8,],col = "0.7")) +geom_line(aes(T,k[8,],col = "0.7")) +
  geom_point(aes(T,k[9,],col = "0.8")) +geom_line(aes(T,k[9,],col = "0.8")) +
  xlab("Time") + 
  ylab("Value of option")+
  ggtitle("Lambda2")


Q2 = c (3412.33,3685.29,3924.87,4109.57,4276,4429.96,3573.7,3873.56,4143.92,4326.02,4526.75,4707.44,3691.5,4063.45,4342.43,4576.15,4808.69,4948.31,3844.65,4246.74,4555.96,4807.54,4989.89,5206.71,3979.05,4362.9,4733.53,4989.49,5211.04,5442.81,4123.71,4551.2,4968.47,5227.06,5410.38,5640.32,4289.68,4717.52,5069.93,5435.65,5647.61,5852.12,4425.51,4910.07,5319.47,5601.85,5833.69,6052.36,4552.49,5063.87,5456.65,5750.44,6018.35,6228.91)
lambda1 = c(seq(0,0.4 ,by = 0.08))

k = matrix(0, nrow = 9, ncol = 6)

for (i in seq(1,54,6)){
  print(i)
  k[(i-1)/6 +1,] = Q2[i:(i+5)]
}

ggplot() + geom_point(aes(T,k[1,],col = "0")) +geom_line(aes(T,k[1,],col = "0")) +
  geom_point(aes(T,k[2,],col = "0.05")) +geom_line(aes(T,k[2,],col = "0.05")) +
  geom_point(aes(T,k[3,],col = "0.1")) +geom_line(aes(T,k[3,],col = "0.1")) +
  geom_point(aes(T,k[4,],col = "0.15")) +geom_line(aes(T,k[4,],col = "0.15")) +
  geom_point(aes(T,k[5,],col = "0.2")) +geom_line(aes(T,k[5,],col = "0.2")) +
  geom_point(aes(T,k[6,],col = "0.25")) +geom_line(aes(T,k[6,],col = "0.25")) +
  geom_point(aes(T,k[7,],col = "0.3")) +geom_line(aes(T,k[7,],col = "0.3")) +
  geom_point(aes(T,k[8,],col = "0.35")) +geom_line(aes(T,k[8,],col = "0.35")) +
  geom_point(aes(T,k[9,],col = "0.4")) +geom_line(aes(T,k[9,],col = "0.4")) +
  xlab("T") + 
  ylab("Value of option")+
  ggtitle("Lambda1")



Q1 = c (0.2957,0.4404,0.5725,0.6746,0.7635,0.8288,0.4854,0.6394,0.7556,0.8372,0.889,0.9296,0.6253,0.7547,0.8608,0.9172,0.954,0.9727,0.727,0.849,0.9238,0.9557,0.9766,0.9872,0.8022,0.9004,0.955,0.9764,0.9901,0.9959,0.8595,0.9333,0.9742,0.9898,0.9959,0.9987,0.8981,0.9574,0.9865,0.9952,0.9979,0.9992,0.9209,0.9741,0.9909,0.9978,0.9996,0.9996,0.9431,0.9851,0.9959,0.9989,0.9994,0.9998)
T = c(seq(3,8 ,by = 1))

k = matrix(0, nrow = 9, ncol = 6)

for (i in seq(1,54,6)){
  print(i)
  k[(i-1)/6 +1,] = Q1[i:(i+5)]
}

ggplot() + geom_point(aes(T,k[1,],col = "0")) +geom_line(aes(T,k[1,],col = "0")) +
  geom_point(aes(T,k[2,],col = "0.1")) +geom_line(aes(T,k[2,],col = "0.1")) +
  geom_point(aes(T,k[3,],col = "0.2")) +geom_line(aes(T,k[3,],col = "0.2")) +
  geom_point(aes(T,k[4,],col = "0.3")) +geom_line(aes(T,k[4,],col = "0.3")) +
  geom_point(aes(T,k[5,],col = "0.4")) +geom_line(aes(T,k[5,],col = "0.4")) +
  geom_point(aes(T,k[6,],col = "0.5")) +geom_line(aes(T,k[6,],col = "0.5")) +
  geom_point(aes(T,k[7,],col = "0.6")) +geom_line(aes(T,k[7,],col = "0.6")) +
  geom_point(aes(T,k[8,],col = "0.7")) +geom_line(aes(T,k[8,],col = "0.7")) +
  geom_point(aes(T,k[9,],col = "0.8")) +geom_line(aes(T,k[9,],col = "0.8")) +
  xlab("Time") + 
  ylab("Probability of Default")+
  ggtitle("Lambda2")




Q1 = c (0.7326,0.8474,0.9195,0.9564,0.9791,0.9891,0.7507,0.8628,0.929,0.9627,0.9827,0.9913,0.7687,0.8762,0.9382,0.9694,0.9862,0.9935,0.7894,0.8906,0.9482,0.9737,0.987,0.9936,0.8013,0.9004,0.955,0.9764,0.9901,0.9959,0.8187,0.9112,0.9615,0.9795,0.9913,0.9966,0.8342,0.9159,0.965,0.985,0.9944,0.9977,0.8432,0.9298,0.9685,0.988,0.9946,0.9988,0.8585,0.9341,0.9723,0.99,0.9951,0.9982)
T = c(seq(3,8 ,by = 1))

k = matrix(0, nrow = 9, ncol = 6)

for (i in seq(1,54,6)){
  print(i)
  k[(i-1)/6 +1,] = Q1[i:(i+5)]
}

ggplot() + geom_point(aes(T,k[1,],col = "0")) +geom_line(aes(T,k[1,],col = "0")) +
  geom_point(aes(T,k[2,],col = "0.05")) +geom_line(aes(T,k[2,],col = "0.05")) +
  geom_point(aes(T,k[3,],col = "0.1")) +geom_line(aes(T,k[3,],col = "0.1")) +
  geom_point(aes(T,k[4,],col = "0.15")) +geom_line(aes(T,k[4,],col = "0.15")) +
  geom_point(aes(T,k[5,],col = "0.2")) +geom_line(aes(T,k[5,],col = "0.2")) +
  geom_point(aes(T,k[6,],col = "0.25")) +geom_line(aes(T,k[6,],col = "0.25")) +
  geom_point(aes(T,k[7,],col = "0.3")) +geom_line(aes(T,k[7,],col = "0.3")) +
  geom_point(aes(T,k[8,],col = "0.35")) +geom_line(aes(T,k[8,],col = "0.35")) +
  geom_point(aes(T,k[9,],col = "0.4")) +geom_line(aes(T,k[9,],col = "0.4")) +
  xlab("Time") + 
  ylab("Probability of Default")+
  ggtitle("Lambda1")




Q1 = c (14.3566,14.3224,13.9329,13.0051,12.4436,11.6971,19.296,17.4742,15.5841,14.1789,12.8183,11.5066,19.8461,17.5598,15.2944,13.4664,11.82,10.4548,19.6996,17.1827,14.7402,12.2198,10.6795,9.39219,19.1535,16.1954,13.3823,11.265,9.37762,8.08297,18.2998,14.9428,12.1979,10.0174,8.53063,7.32653,17.1024,14.0722,11.2157,9.20558,7.67491,6.58155,16.2845,12.9209,10.2697,8.4166,6.9746,6.02055,15.6462,11.8502,9.47276,7.62999,6.43207,5.4938)
T = c(seq(3,8 ,by = 1))

k = matrix(0, nrow = 9, ncol = 6)

for (i in seq(1,54,6)){
  print(i)
  k[(i-1)/6 +1,] = Q1[i:(i+5)]
}

ggplot() + geom_point(aes(T,k[1,],col = "0")) +geom_line(aes(T,k[1,],col = "0")) +
  geom_point(aes(T,k[2,],col = "0.1")) +geom_line(aes(T,k[2,],col = "0.1")) +
  geom_point(aes(T,k[3,],col = "0.2")) +geom_line(aes(T,k[3,],col = "0.2")) +
  geom_point(aes(T,k[4,],col = "0.3")) +geom_line(aes(T,k[4,],col = "0.3")) +
  geom_point(aes(T,k[5,],col = "0.4")) +geom_line(aes(T,k[5,],col = "0.4")) +
  geom_point(aes(T,k[6,],col = "0.5")) +geom_line(aes(T,k[6,],col = "0.5")) +
  geom_point(aes(T,k[7,],col = "0.6")) +geom_line(aes(T,k[7,],col = "0.6")) +
  geom_point(aes(T,k[8,],col = "0.7")) +geom_line(aes(T,k[8,],col = "0.7")) +
  geom_point(aes(T,k[9,],col = "0.8")) +geom_line(aes(T,k[9,],col = "0.8")) +
  xlab("Time") + 
  ylab("Probability of Default")+
  ggtitle("Lambda2")




Q1 = c (22.3529,18.9868,16.1287,13.7951,11.5653,10.0156,21.4,18.1445,15.302,13.1184,10.939,9.47835,20.5063,17.1646,14.5246,12.4097,10.3645,8.9277,19.8238,16.7414,13.9809,11.8054,9.89213,8.4424,19.1066,16.1954,13.3823,11.265,9.37762,8.08297,18.2843,15.3833,12.7851,10.7353,9.08272,7.72426,17.8647,14.7354,12.411,10.281,8.64244,7.4304,17.0552,14.2536,11.8405,9.92907,8.31189,7.06784,16.7467,13.8028,11.3822,9.48972,7.86131,6.71916)
T = c(seq(3,8 ,by = 1))

k = matrix(0, nrow = 9, ncol = 6)

for (i in seq(1,54,6)){
  print(i)
  k[(i-1)/6 +1,] = Q1[i:(i+5)]
}

ggplot() + geom_point(aes(T,k[1,],col = "0")) +geom_line(aes(T,k[1,],col = "0")) +
  geom_point(aes(T,k[2,],col = "0.05")) +geom_line(aes(T,k[2,],col = "0.05")) +
  geom_point(aes(T,k[3,],col = "0.1")) +geom_line(aes(T,k[3,],col = "0.1")) +
  geom_point(aes(T,k[4,],col = "0.15")) +geom_line(aes(T,k[4,],col = "0.15")) +
  geom_point(aes(T,k[5,],col = "0.2")) +geom_line(aes(T,k[5,],col = "0.2")) +
  geom_point(aes(T,k[6,],col = "0.25")) +geom_line(aes(T,k[6,],col = "0.25")) +
  geom_point(aes(T,k[7,],col = "0.3")) +geom_line(aes(T,k[7,],col = "0.3")) +
  geom_point(aes(T,k[8,],col = "0.35")) +geom_line(aes(T,k[8,],col = "0.35")) +
  geom_point(aes(T,k[9,],col = "0.4")) +geom_line(aes(T,k[9,],col = "0.4")) +
  xlab("Time") + 
  ylab("Expected Default Time")+
  ggtitle("Lambda1")


call = c(9.26829,12.557,15.9365,19.4012,22.9483,26.5772,30.2871,34.0784,37.952,41.9082)
vol = c(seq(0.12, 0.48,0.04))
ggplot() + geom_point(aes(vol,call,col = "call")) +geom_line(aes(vol,call,col = "call")) +
  xlab("Volitility") + 
  ylab("call option price")

put = c(9.49876,12.2981,15.0386,17.7167,20.3307,22.8816,25.3698,27.7957,30.1596,32.4626)
vol = c(seq(0.12, 0.48,0.04))
ggplot() + geom_point(aes(vol,put,col = "put")) +geom_line(aes(vol,col = "put")) +
  xlab("Volitility") + 
  ylab("call option price")
