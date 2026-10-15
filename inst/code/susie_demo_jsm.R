# Script to test out the new SuSiE interface using a data set that I
# generated for my JSM tutorial.
library(ebverse)
set.seed(1)
load("../datasets/jsm_data.RData")
fit <- fit_susie(X,y,L = 2,numiter = 10,verbose = TRUE)
plot(fit$PIP[,1],pch = 20,col = "dodgerblue",ylim = c(0,1),
     xlab = "SNP",ylab = "PIP")
points(fit$PIP[,2],pch = 20,col = "red")
