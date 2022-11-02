library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)
library(stringr)

#gname = c("theta_chain.eps",sep="")  
postscript(gname_theta,width=12,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,3),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(Theta_chain[1,], type = 'l', xlab = '', ylab = '', lwd = 1)
#plot(Theta_chain[2,], type = 'l', xlab = '', ylab = '', lwd = 1)
#plot(Theta_chain[3,], type = 'l', xlab = '', ylab = '', lwd = 1)
#plot(Theta_chain[4,], type = 'l', xlab = '', ylab = '', lwd = 1)
mtext(expression(Theta[11]), side = 1, line= 2)
plot(Theta_chain[14,], type = 'l', xlab = '', ylab = '', lwd = 1)
#plot(Theta_chain[15,], type = 'l', xlab = '', ylab = '', lwd = 1)
#plot(Theta_chain[16,], type = 'l', xlab = '', ylab = '', lwd = 1)
#plot(Theta_chain[17,], type = 'l', xlab = '', ylab = '', lwd = 1)
mtext(expression(Theta[21]), side = 1, line= 2)
plot(Theta_chain[27,], type = 'l', xlab = '', ylab = '', lwd = 1)
#plot(Theta_chain[28,], type = 'l', xlab = '', ylab = '', lwd = 1)
#plot(Theta_chain[29,], type = 'l', xlab = '', ylab = '', lwd = 1)
#plot(Theta_chain[30,], type = 'l', xlab = '', ylab = '', lwd = 1)
mtext(expression(Theta[31]), side = 1, line= 2)
dev.off()

#gname = c("xi_chain.eps",sep="")  
postscript(gname_xi,width=12,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,3),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(xi_chain[1,], type = 'l', xlab = '', ylab = '', lwd = 1)
abline(h=xi[1], lty = 2, col = 2, lwd = 2)
mtext(expression(xi[11]), side = 1, line= 2)
plot(xi_chain[11,], type = 'l', xlab = '', ylab = '', lwd = 1)
abline(h=xi[11], lty = 2, col = 2, lwd = 2)
mtext(expression(xi[42]), side = 1, line= 2)
plot(xi_chain[27,], type = 'l', xlab = '', ylab = '', lwd = 1)
abline(h=xi[27], lty = 2, col = 2, lwd = 2)
mtext(expression(xi[93]), side = 1, line= 2)
#plot(xi_chain[41,], type = 'l', xlab = '', ylab = '', lwd = 1)
#abline(h=xi[41], lty = 2, col = 2, lwd = 2)
dev.off()

#gname = c("tau_chain.eps",sep="")  
postscript(gname_tau,width=12,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,3),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(tau_chain[1,], type = 'l', xlab = '', ylab = '', lwd = 1)
abline(h=tau[1], lty = 2, col = 2, lwd = 2)
mtext(expression(tau[1]), side = 1, line= 2)
plot(tau_chain[2,], type = 'l', xlab = '', ylab = '', lwd = 1)
abline(h=d[2], lty = 2, col = 2, lwd = 2)
mtext(expression(d[2]), side = 1, line= 2)
plot(tau_chain[3,], type = 'l', xlab = '', ylab = '', lwd = 1)
abline(h=d[3], lty = 2, col = 2, lwd = 2)
mtext(expression(d[3]), side = 1, line= 2)
dev.off()

#gname = c("Lambda_chain.eps",sep="")  
postscript(gname_Lambda,width=12,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,3),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(sqrt(Lambda_chain[1,]), type = 'l', xlab = '', ylab = '', lwd = 1)
abline(h=lambda[1], lty = 2, col = 2, lwd = 2)
mtext(expression(lambda[1]), side = 1, line= 2)
plot(sqrt(Lambda_chain[2,]), type = 'l', xlab = '', ylab = '', lwd = 1)
abline(h=lambda[2], lty = 2, col = 2, lwd = 2)
mtext(expression(lambda[2]), side = 1, line= 2)
plot(sqrt(Lambda_chain[3,]), type = 'l', xlab = '', ylab = '', lwd = 1)
abline(h=lambda[3], lty = 2, col = 2, lwd = 2)
mtext(expression(lambda[3]), side = 1, line= 2)
dev.off()

#gname = c("sig_chain.eps",sep="")  
#postscript(gname,width=6,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
#par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

#plot(sig2e_chain, type = 'l', xlab = '', ylab = '', lwd = 1)
#mtext(expression(sigma^2), side = 1, line= 2)

#dev.off()

#Lambda_mean <- apply(Lambda_chain[,-(1:2000)], 1, mean)
#Theta_mean <- matrix(apply(Theta_chain[ ,-(1:2000)], 1, mean), nc = K)
#mat <- Theta_mean%*%diag(Lambda_mean)%*%t(Theta_mean)
#Theta_mean <- eigen(mat)$vectors[,1:K]


#gname = c("phi_est.eps",sep="")  
postscript(gname_phi,width=12,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,3),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(phi_est[1,], type = 'l', ylim = c(-1.5,1.5), xlab = '', ylab = '', lwd = 2)
lines(phi[1,], lty = 2, col = 2, lwd = 2)
mtext(expression(phi[1](t)), side = 1, line= 2)

plot(phi_est[2,], type = 'l', ylim = c(-1.5,1.5), xlab = '', ylab = '', lwd = 2)
lines(phi[2,], lty = 2, col = 2, lwd = 2)
mtext(expression(phi[2](t)), side = 1, line= 2)

plot(phi_est[3,], type = 'l', ylim = c(-1.5,1.5), xlab = '', ylab = '', lwd = 2)
lines(phi[3,], lty = 2, col = 2, lwd = 2)
mtext(expression(phi[3](t)), side = 1, line= 2)
dev.off()


#est_curve <- basismat%*%apply(c_chain[,-(1:2000)], 1, mean)

#gname = c("mean_est.eps",sep="")  
postscript(gname_mean,width=5,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(est_curve, type = 'l', ylab = '', xlab = '', lwd = 2)
lines(meancurve, lty = 2, col = 2, lwd = 2)
mtext(expression(mu(t)), side = 1, line= 2)
dev.off()


