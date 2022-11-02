
#gname = c("theta_chain.eps",sep="")  
postscript(gname_theta,width=12,height=15,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(6,2),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(Theta_chain[1,], type = 'l', xlab = '', ylab = '', lwd = 1)
plot(Theta_chain[2,], type = 'l', xlab = '', ylab = '', lwd = 1)
plot(Theta_chain[3,], type = 'l', xlab = '', ylab = '', lwd = 1)
plot(Theta_chain[4,], type = 'l', xlab = '', ylab = '', lwd = 1)

plot(Theta_chain[14,], type = 'l', xlab = '', ylab = '', lwd = 1)
plot(Theta_chain[15,], type = 'l', xlab = '', ylab = '', lwd = 1)
plot(Theta_chain[16,], type = 'l', xlab = '', ylab = '', lwd = 1)
plot(Theta_chain[17,], type = 'l', xlab = '', ylab = '', lwd = 1)

plot(Theta_chain[27,], type = 'l', xlab = '', ylab = '', lwd = 1)
plot(Theta_chain[28,], type = 'l', xlab = '', ylab = '', lwd = 1)
plot(Theta_chain[29,], type = 'l', xlab = '', ylab = '', lwd = 1)
plot(Theta_chain[30,], type = 'l', xlab = '', ylab = '', lwd = 1)

dev.off()

#gname = c("xi_chain.eps",sep="")  
postscript(gname_xi,width=12,height=5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(2,2),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(xi_chain[1,], type = 'l', xlab = '', ylab = '', lwd = 1)
abline(h=xi[1], lty = 2, col = 2, lwd = 2)
plot(xi_chain[11,], type = 'l', xlab = '', ylab = '', lwd = 1)
abline(h=xi[11], lty = 2, col = 2, lwd = 2)
plot(xi_chain[31,], type = 'l', xlab = '', ylab = '', lwd = 1)
abline(h=xi[31], lty = 2, col = 2, lwd = 2)
plot(xi_chain[41,], type = 'l', xlab = '', ylab = '', lwd = 1)
abline(h=xi[41], lty = 2, col = 2, lwd = 2)
dev.off()

#gname = c("Lambda_chain.eps",sep="")  
postscript(gname_Lambda,width=12,height=5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(2,2),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

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


#gname = c("phi_est.eps",sep="")  
postscript(gname_phi,width=8,height=5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(2,2),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

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


#gname = c("mean_est.eps",sep="")  
postscript(gname_mean,width=5,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(est_curve, type = 'l', ylab = '', xlab = '', lwd = 2)
lines(meancurve, lty = 2, col = 2, lwd = 2)
mtext(expression(mu(t)), side = 1, line= 2)
dev.off()

#gname = c("tau_chain.eps",sep="")  
postscript(gname_tau,width=12,height=5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(2,2),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(tau_chain[1,], type = 'l', xlab = '', ylab = '', lwd = 1)
plot(tau_chain[2,], type = 'l', xlab = '', ylab = '', lwd = 1)
plot(tau_chain[3,], type = 'l', xlab = '', ylab = '', lwd = 1)
dev.off()


