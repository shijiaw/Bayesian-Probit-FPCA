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


