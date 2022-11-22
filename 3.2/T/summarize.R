library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)

nrep <- 30
K <- 3
T21 <- data.frame(T = rep(21, nrep))
T21$mu = unlist(read.table("T21/meancurve_RMSERes.txt", header = FALSE))
phi = matrix(unlist(read.table("T21/phi_RMSERes.txt", header = FALSE)), nc = K)
T21$phi1 <- phi[,1]
T21$phi2 <- phi[,2]
T21$phi3 <- phi[,3]
lambda = matrix(unlist(read.table("T21/Lambda_mean.txt", header = FALSE)), nc = K)
T21$lambda1 <- lambda[,1]
T21$lambda2 <- lambda[,2]
T21$lambda3 <- lambda[,3]
tau = matrix(unlist(read.table("T21/tau_mean.txt", header = FALSE)), nc = K)
T21$tau1 <- tau[,1]
T21$tau2 <- tau[,2]
T21$tau3 <- tau[,3]

T101 <- data.frame(T = rep(101, nrep))
T101$mu = unlist(read.table("T101/meancurve_RMSERes.txt", header = FALSE))
phi = matrix(unlist(read.table("T101/phi_RMSERes.txt", header = FALSE)), nc = K)
T101$phi1 <- phi[,1]
T101$phi2 <- phi[,2]
T101$phi3 <- phi[,3]
lambda = matrix(unlist(read.table("T101/Lambda_mean.txt", header = FALSE)), nc = K)
T101$lambda1 <- lambda[,1]
T101$lambda2 <- lambda[,2]
T101$lambda3 <- lambda[,3]
tau = matrix(unlist(read.table("T101/tau_mean.txt", header = FALSE)), nc = K)
T101$tau1 <- tau[,1]
T101$tau2 <- tau[,2]
T101$tau3 <- tau[,3]

Results <- rbind(T21, T101)
Results$T = as.factor(Results$T)

p0 <- ggplot(Results, aes(T, mu))
p1 <- ggplot(Results, aes(T, phi1))
p2 <- ggplot(Results, aes(T, phi2))
p3 <- ggplot(Results, aes(T, phi3))
p4 <- ggplot(Results, aes(T, lambda1))
p5 <- ggplot(Results, aes(T, lambda2))
p6 <- ggplot(Results, aes(T, lambda3))
p7 <- ggplot(Results, aes(T, tau1))
p8 <- ggplot(Results, aes(T, tau2))
p9 <- ggplot(Results, aes(T, tau3))

gname = c("T.eps",sep="")  
postscript(gname,width=10,height=5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

ggarrange(p0 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = T))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(mu))),
          p1 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = T))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(phi[1]))),
          p2 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = T))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(phi[2]))),
          p3 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = T))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(phi[3]))),
          p4 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = T))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(lambda[1]))+
            geom_hline(yintercept=3, linetype="dotted"),
          p5 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = T))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(lambda[2]))+
            geom_hline(yintercept=1.5, linetype="dotted"),
          p6 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = T))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(lambda[3]))+
            geom_hline(yintercept=0.75, linetype="dotted"),
          p7 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = T))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(tau[1]))+
            geom_hline(yintercept=-3.5, linetype="dotted"),
          p8 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = T))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(gamma[2]))+
            geom_hline(yintercept=0.5, linetype="dotted"),
          p9 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = T))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(gamma[3]))+
            geom_hline(yintercept=1, linetype="dotted"),
          #p4 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = sigma))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab("RMSE"),
          ncol = 5, nrow = 2, common.legend = TRUE
          )

dev.off()


