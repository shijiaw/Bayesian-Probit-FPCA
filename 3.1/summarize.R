library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)

nrep <- 30
K <- 3

probit <- unlist(read.table("model1.txt", header = FALSE))
efpca <- unlist(read.table("model2.txt", header = FALSE))
face <- unlist(read.table("model3.txt", header = FALSE))
sc <- unlist(read.table("model4.txt", header = FALSE))
percent <- data.frame(result = c(probit, efpca, face, sc), method = rep(c('probit', 'efpca', 'face', 'sc'), each = nrep))
percent$method = as.factor(percent$method)

p0 <- ggplot(percent, aes(method, result))

gname = c("S1.eps",sep="")  
postscript(gname,width=3.5,height=3.5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

ggarrange(p0 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = method))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab("Classification rate"),
          ncol = 1, nrow = 1, common.legend = TRUE
          )

dev.off()


probit <- unlist(read.table("meancurve_RMSE.txt", header = FALSE))
efpca <- unlist(read.table("fpcamean_RMSE.txt", header = FALSE))
face <- unlist(read.table("facemean_RMSE.txt", header = FALSE))
sc <- unlist(read.table("scmean_RMSE.txt", header = FALSE))

meancurve <- data.frame(RMSE = c(probit, efpca, face, sc), method = rep(c('probit', 'efpca', 'face', 'sc'), each = nrep))
meancurve$method = as.factor(meancurve$method)


gname = c("S2.eps",sep="")  
postscript(gname,width=3.5,height=3.5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

ggarrange(p1 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = method))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab("Method"),
          ncol = 1, nrow = 1, common.legend = TRUE
)

dev.off()


probit_fpc <- matrix(unlist(read.table("phi_RMSE.txt", header = FALSE)), ncol = K)
efpca_fpc <- matrix(unlist(read.table("fpcafpc_RMSE.txt", header = FALSE)), ncol = K)
face_fpc <- matrix(unlist(read.table("facefpc_RMSE.txt", header = FALSE)), ncol = K)
sc_fpc <- matrix(unlist(read.table("scfpc_RMSE.txt", header = FALSE)), ncol = K)

fpc_RMSE1 <- data.frame(RMSE = c(probit_fpc[,1], efpca_fpc[,1], face_fpc[,1], sc_fpc[,1]), method = rep(c('probit', 'efpca', 'face', 'sc'), each = nrep))
fpc_RMSE1$method = as.factor(fpc_RMSE1$method)

fpc_RMSE2 <- data.frame(RMSE = c(probit_fpc[,2], efpca_fpc[,2], face_fpc[,2], sc_fpc[,2]), method = rep(c('probit', 'efpca', 'face', 'sc'), each = nrep))
fpc_RMSE2$method = as.factor(fpc_RMSE2$method)

fpc_RMSE3 <- data.frame(RMSE = c(probit_fpc[,3], efpca_fpc[,3], face_fpc[,3], sc_fpc[,3]), method = rep(c('probit', 'efpca', 'face', 'sc'), each = nrep))
fpc_RMSE3$method = as.factor(fpc_RMSE3$method)

p1 <- ggplot(meancurve, aes(method, RMSE))+ geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = method))+ theme_bw()+rremove("x.text") +xlab(expression(mu(t)))
p2 <- ggplot(fpc_RMSE1, aes(method, RMSE))+ geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = method))+ theme_bw()+rremove("x.text")+rremove("ylab")+xlab(expression(phi[1](t)))
p3 <- ggplot(fpc_RMSE2, aes(method, RMSE))+ geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = method))+ theme_bw()+rremove("x.text")+rremove("ylab")+xlab("mean curve")+xlab(expression(phi[2](t)))
p4 <- ggplot(fpc_RMSE3, aes(method, RMSE))+ geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = method))+ theme_bw()+rremove("x.text")+rremove("ylab")+xlab("mean curve")+xlab(expression(phi[3](t)))

gname = c("S2.eps",sep="")  
postscript(gname,width=10,height=3.5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
ggarrange(
  p1,
  p2,
  p3,
  p4,
  ncol = 4, nrow = 1, common.legend = TRUE, labels = NULL
)
dev.off()


