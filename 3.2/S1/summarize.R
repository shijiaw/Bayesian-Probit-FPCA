library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)

nrep <- 30

dic_temp <- matrix(as.numeric(read.table("dic_Res.txt", sep = ' ')), nc = 3)
epd_temp <- matrix(as.numeric(read.table("epd_Res.txt", sep = ' ')), nc = 3)
eaic_temp <- matrix(as.numeric(read.table("eaic_Res.txt", sep = ' ')), nc = 3)
waic_temp <- matrix(as.numeric(read.table("waic_Res.txt", sep = ' ')), nc = 3)
ebic_temp <- matrix(as.numeric(read.table("ebic_Res.txt", sep = ' ')), nc = 3)


dic <- as.vector(apply(dic_temp, 1, function(x)(x- x[2])))
epd <- as.vector(apply(epd_temp, 1, function(x)(x- x[2])))
eaic <- as.vector(apply(eaic_temp, 1, function(x)(x- x[2])))
waic <- as.vector(apply(waic_temp, 1, function(x)(x- x[2])))
ebic <- as.vector(apply(ebic_temp, 1, function(x)(x- x[2])))

S1 <- data.frame(K = rep(2:4, 30))
S1$dic <- dic
S1$epd <- epd
S1$eaic <- eaic
S1$ebic <- ebic
S1$waic <- waic
S1$K <- as.factor(S1$K)
p0 <- ggplot(S1, aes(K, dic)) 
p1 <- ggplot(S1, aes(K, waic)) 
p2 <- ggplot(S1, aes(K, eaic)) 
p3 <- ggplot(S1, aes(K, ebic)) 
p4 <- ggplot(S1, aes(K, epd)) 


gname = c("dicS1.eps",sep="")  
postscript(gname,width=10,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

ggarrange(p0 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = K))+ theme_bw()+rremove("y.text")+rremove("xlab")+ ylab(expression(DIC-DIC[3]))+
            geom_hline(yintercept=0, linetype="dotted"),
          p1 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = K))+ theme_bw()+rremove("y.text")+rremove("xlab")+ ylab(expression(WAIC-WAIC[3]))+
            geom_hline(yintercept=0, linetype="dotted"),
          p2 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = K))+ theme_bw()+rremove("y.text")+rremove("xlab")+ ylab(expression(EAIC-EAIC[3]))+
            geom_hline(yintercept=0, linetype="dotted"),
          p3 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = K))+ theme_bw()+rremove("y.text")+rremove("xlab")+ ylab(expression(EBIC-EBIC[3]))+
            geom_hline(yintercept=0, linetype="dotted"),
          p4 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = K))+ theme_bw()+rremove("y.text")+rremove("xlab")+ ylab(expression(EPD-EPD[3]))+
            geom_hline(yintercept=0, linetype="dotted"),
          ncol = 5, nrow = 1, common.legend = TRUE
          )

dev.off()


