library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)

nrep <- 30
K <- 3
N50 <- data.frame(N = rep(100, nrep))
N50$mu = unlist(read.table("N50/meancurve_RMSERes.txt", header = FALSE))
phi = matrix(unlist(read.table("N50/phi_RMSERes.txt", header = FALSE)), nc = K)
N50$phi1 <- phi[,1]
N50$phi2 <- phi[,2]
N50$phi3 <- phi[,3]
lambda = matrix(unlist(read.table("N50/Lambda_mean.txt", header = FALSE)), nc = K)
N50$lambda1 <- lambda[,1]
N50$lambda2 <- lambda[,2]
N50$lambda3 <- lambda[,3]
tau = matrix(unlist(read.table("N50/tau_mean.txt", header = FALSE)), nc = K)
N50$tau1 <- tau[,1]
N50$tau2 <- tau[,2]
N50$tau3 <- tau[,3]

N500 <- data.frame(N = rep(500, nrep))
N500$mu = unlist(read.table("N500/meancurve_RMSERes.txt", header = FALSE))
phi = matrix(unlist(read.table("N500/phi_RMSERes.txt", header = FALSE)), nc = K)
N500$phi1 <- phi[,1]
N500$phi2 <- phi[,2]
N500$phi3 <- phi[,3]
lambda = matrix(unlist(read.table("N500/Lambda_mean.txt", header = FALSE)), nc = K)
N500$lambda1 <- lambda[,1]
N500$lambda2 <- lambda[,2]
N500$lambda3 <- lambda[,3]
tau = matrix(unlist(read.table("N500/tau_mean.txt", header = FALSE)), nc = K)
N500$tau1 <- tau[,1]
N500$tau2 <- tau[,2]
N500$tau3 <- tau[,3]

N2000 <- data.frame(N = rep(2000, nrep))
N2000$mu = unlist(read.table("N2000/meancurve_RMSERes.txt", header = FALSE))
phi = matrix(unlist(read.table("N2000/phi_RMSERes.txt", header = FALSE)), nc = K)
N2000$phi1 <- phi[,1]
N2000$phi2 <- phi[,2]
N2000$phi3 <- phi[,3]
lambda = matrix(unlist(read.table("N2000/Lambda_mean.txt", header = FALSE)), nc = K)
N2000$lambda1 <- lambda[,1]
N2000$lambda2 <- lambda[,2]
N2000$lambda3 <- lambda[,3]
tau = matrix(unlist(read.table("N2000/tau_mean.txt", header = FALSE)), nc = K)
N2000$tau1 <- tau[,1]
N2000$tau2 <- tau[,2]
N2000$tau3 <- tau[,3]

Results <- rbind(N50, N500, N2000)
Results$N = as.factor(Results$N)

#Resultssummary <- summarySE(Results, measurevar=c("nu"), groupvars=c("sigma"))
p0 <- ggplot(Results, aes(N, mu))
p1 <- ggplot(Results, aes(N, phi1))
p2 <- ggplot(Results, aes(N, phi2))
p3 <- ggplot(Results, aes(N, phi3))
p4 <- ggplot(Results, aes(N, lambda1))
p5 <- ggplot(Results, aes(N, lambda2))
p6 <- ggplot(Results, aes(N, lambda3))
p7 <- ggplot(Results, aes(N, tau1))
p8 <- ggplot(Results, aes(N, tau2))
p9 <- ggplot(Results, aes(N, tau3))

gname = c("N.eps",sep="")  
postscript(gname,width=10,height=5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

ggarrange(p0 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = N))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(mu))),
          p1 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = N))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(phi[1]))),
          p2 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = N))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(phi[2]))),
          p3 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = N))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(phi[3]))),
          p4 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = N))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(lambda[1]))+
            geom_hline(yintercept=3, linetype="dotted"),
          p5 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = N))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(lambda[2]))+
            geom_hline(yintercept=1.5, linetype="dotted"),
          p6 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = N))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(lambda[3]))+
            geom_hline(yintercept=0.75, linetype="dotted"),
          p7 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = N))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(tau[1]))+
            geom_hline(yintercept=-3.5, linetype="dotted"),
          p8 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = N))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(gamma[2]))+
            geom_hline(yintercept=0.5, linetype="dotted"),
          p9 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = N))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(gamma[3]))+
            geom_hline(yintercept=1, linetype="dotted"),
          #p4 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = sigma))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab("RMSE"),
          ncol = 5, nrow = 2, common.legend = TRUE
          )

dev.off()


