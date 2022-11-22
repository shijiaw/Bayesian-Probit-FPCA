library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)

nrep <- 30
K <- 3
tau1 <- data.frame(threshold = rep(1, nrep))
tau1$mu = unlist(read.table("tau1/meancurve_RMSERes.txt", header = FALSE))
phi = matrix(unlist(read.table("tau1/phi_RMSERes.txt", header = FALSE)), nc = K)
tau1$phi1 <- phi[,1]
tau1$phi2 <- phi[,2]
tau1$phi3 <- phi[,3]
lambda = matrix(unlist(read.table("tau1/Lambda_mean.txt", header = FALSE)), nc = K)
tau1$lambda1 <- lambda[,1]
tau1$lambda2 <- lambda[,2]
tau1$lambda3 <- lambda[,3]
tau = matrix(unlist(read.table("tau1/tau_mean.txt", header = FALSE)), nc = K)
tau1$tau1 <- tau[,1] - (- 5)
tau1$tau2 <- tau[,2] - 0.5
tau1$tau3 <- tau[,3] - 2


tau2 <- data.frame(threshold = rep(2, nrep))
tau2$mu = unlist(read.table("tau2/meancurve_RMSERes.txt", header = FALSE))
phi = matrix(unlist(read.table("tau2/phi_RMSERes.txt", header = FALSE)), nc = K)
tau2$phi1 <- phi[,1]
tau2$phi2 <- phi[,2]
tau2$phi3 <- phi[,3]
lambda = matrix(unlist(read.table("tau2/Lambda_mean.txt", header = FALSE)), nc = K)
tau2$lambda1 <- lambda[,1]
tau2$lambda2 <- lambda[,2]
tau2$lambda3 <- lambda[,3]
tau = matrix(unlist(read.table("tau2/tau_mean.txt", header = FALSE)), nc = K)
tau2$tau1 <- tau[,1] - (-5)
tau2$tau2 <- tau[,2] - 1
tau2$tau3 <- tau[,3] - (-1)

tau3 <- data.frame(threshold = rep(3, nrep))
tau3$mu = unlist(read.table("tau3/meancurve_RMSERes.txt", header = FALSE))
phi = matrix(unlist(read.table("tau3/phi_RMSERes.txt", header = FALSE)), nc = K)
tau3$phi1 <- phi[,1]
tau3$phi2 <- phi[,2]
tau3$phi3 <- phi[,3]
lambda = matrix(unlist(read.table("tau3/Lambda_mean.txt", header = FALSE)), nc = K)
tau3$lambda1 <- lambda[,1]
tau3$lambda2 <- lambda[,2]
tau3$lambda3 <- lambda[,3]
tau = matrix(unlist(read.table("tau3/tau_mean.txt", header = FALSE)), nc = K)
tau3$tau1 <- tau[,1] - (-3.5)
tau3$tau2 <- tau[,2] - 0.5
tau3$tau3 <- tau[,3] - 1

Results <- rbind(tau1, tau2, tau3)
Results$threshold = as.factor(Results$threshold)

#Resultssummary <- summarySE(Results, measurevar=c("nu"), groupvars=c("sigma"))
p0 <- ggplot(Results, aes(threshold, mu))
p1 <- ggplot(Results, aes(threshold, phi1))
p2 <- ggplot(Results, aes(threshold, phi2))
p3 <- ggplot(Results, aes(threshold, phi3))
p4 <- ggplot(Results, aes(threshold, lambda1))
p5 <- ggplot(Results, aes(threshold, lambda2))
p6 <- ggplot(Results, aes(threshold, lambda3))
p7 <- ggplot(Results, aes(threshold, tau1))
p8 <- ggplot(Results, aes(threshold, tau2))
p9 <- ggplot(Results, aes(threshold, tau3))



gname = c("tau.eps",sep="")  
postscript(gname,width=10,height=5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

ggarrange(p0 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = threshold))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(mu))),
          p1 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = threshold))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(phi[1]))),
          p2 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = threshold))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(phi[2]))),
          p3 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = threshold))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(phi[3]))),
          p4 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = threshold))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(lambda[1]))+
            geom_hline(yintercept=3, linetype="dotted"),
          p5 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = threshold))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(lambda[2]))+
            geom_hline(yintercept=1.5, linetype="dotted"),
          p6 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = threshold))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(lambda[3]))+
            geom_hline(yintercept=0.75, linetype="dotted"),
          p7 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = threshold))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(tau[1])),
          p8 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = threshold))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(gamma[2])),
          p9 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = threshold))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(gamma[3])),
          #p4 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = sigma))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab("RMSE"),
          ncol = 5, nrow = 2, common.legend = TRUE
          )

dev.off()


