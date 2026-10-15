# Script to test out the new SuSiE interface using a data set that I
# generated for my JSM tutorial.
#
# library(ebverse)
# remotes::install_local("ebverse",upgrade = "never")
set.seed(1)
load("../datasets/jsm_data.RData")
X <- scale(X,center = TRUE,scale = FALSE)
y <- drop(scale(y,center = TRUE,scale = FALSE))
fit <- fit_susie(X,y,L = 2,numiter = 10,verbose = TRUE)
plot(fit$PIP[,1],pch = 20,col = "dodgerblue",ylim = c(0,1),
     xlab = "SNP",ylab = "PIP")
points(fit$PIP[,2],pch = 20,col = "red")
j <- which(b_true != 0)
pip <- rowSums(fit$PIP)
points(j,pip[j],pch = 1,col = "black",cex = 1.25)
