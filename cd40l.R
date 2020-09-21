setwd('D:/work/work table/CD40L')
library(survival)
library(cmprsk)
library(pec)
library(dplyr)
library(DescTools)
library(rms)
library(maxstat)
library(etm)
cd40h <-mutate(tbl_df(read.csv('cd40lh.CSV')),nrme=as.factor(NRMe_2017),cd40=rs3092936a!=3,tmar=as.factor(tmar)
               ,comb4=(cd40*1+1)*(thbd+1),  cr1=CR>0)
cd40h0 <- filter(cd40h, Statine==0)
cd40h1 <- filter(cd40h, Statine==1)
cd40h11 <- filter(cd40h1, code<4)
cd40h12 <- filter(cd40h1, code==4)
cd40hg <- filter(cd40h0, !is.na(tr))
length(cd40hg$rs3092920q)

summary(factor(cd40h$rs3092936a))
chrh1surv <- npsurv(Surv(nrmm, nrme)~ tr,data=cd40h12)
par(mgp=c(1.5, 1, 0))
m <-survplot(chrh1surv, state = 1, n.risk=TRUE, xlim=c(0, 60),ylim=c(0, 1),conf.int=0.1,time.inc=10,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, 
             col=c("BLACK","#E41A1C","blue"),
             lty = c(1,1,1,1),xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=1)
legend('topright', c("High-risk","Low-risk") , 
       lty = c(1,1), col=c("#E41A1C","BLACK"), bty='n', cex=0.9,lwd = 2, y.intersp=2)

chrh1surv <- npsurv(Surv( DauerOS_2017, Ose)~ tr,data=cd40h0)
par(mgp=c(1.5, 1, 0))
m <-survplot(chrh1surv, state = 1, n.risk=TRUE, xlim=c(0, 120),ylim=c(0, 1),conf.int=0.1,time.inc=20,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, 
             col=c("BLACK","#E41A1C","blue"),
             lty = c(1,1,1,1),xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=1)
legend('topright', c("High-risk","Low-risk") , 
       lty = c(1,1), col=c("#E41A1C","BLACK"), bty='n', cex=0.9,lwd = 2, y.intersp=2)

chrh1surv <- npsurv(Surv(nrmm, nrme)~ tr,data=cd40h11)
par(mgp=c(1.5, 1, 0))
m <-survplot(chrh1surv, state = 1, n.risk=TRUE, xlim=c(0, 120),ylim=c(0, 1),conf.int=0.1,time.inc=20,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, 
             col=c("#377EB8","BLACK","#E41A1C","#A65628"),
             lty = c(1,1,1,1),xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=1)
legend('topright', c("Low-risk","THBD","CD40L",'Both') , 
       lty = c(1,1), col=c("BLACK","#377EB8","#A65628","#E41A1C"), bty='n', cex=0.9,lwd = 2, y.intersp=2)
########################CD40B######################

cd40b <-mutate(tbl_df(read.csv('CD40B.CSV')),nrme=as.factor((CR==2)*1),cd40=rs3092936d,osmm=osm/30.5,
               comb4=(cd40*1+1)*(thbd+1),cr1= CR>0)
cd40bg <- filter(cd40b, !is.na(tr))
cd40b$nrme
chrh1surv <- npsurv(Surv(nrmm, cr1)~ tr,data=cd40bg)
par(mgp=c(1.5, 1, 0))
m <-survplot(chrh1surv, state =1, n.risk=TRUE, xlim=c(0, 120),ylim=c(0, 1),conf.int=0.1,time.inc=20,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, 
             col=c("#377EB8","BLACK","#E41A1C","#A65628"),
             lty = c(1,1,1,1),xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=1)
legend('topright', c("Low-risk","THBD","CD40L",'Both') ,
       lty = c(1,1), col=c("BLACK","#377EB8","#A65628","#E41A1C"), bty='n', cex=0.9,lwd = 2, y.intersp=2)


chrh1surv <- npsurv(Surv(nrmm, cr1)~ tr,data=cd40bg)
par(mgp=c(1.5, 1, 0))
m <-survplot(chrh1surv, state =1, n.risk=TRUE, xlim=c(0, 120),ylim=c(0, 1),conf.int=0.1,time.inc=20,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, 
             col=c("BLACK","#E41A1C"),
             lty = c(1,1,1,1),xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=1)
legend('topright', c("High-risk","Low-risk") , 
       lty = c(1,1), col=c("#E41A1C","BLACK"), bty='n', cex=0.9,lwd = 2, y.intersp=2)

chrh1surv <- npsurv(Surv(osmm, ose)~ tr,data=cd40bg)
par(mgp=c(1.5, 1, 0))
m <-survplot(chrh1surv,  n.risk=TRUE, xlim=c(0, 120),ylim=c(0, 1),conf.int=0.1,time.inc=20,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, 
             col=c("BLACK","#E41A1C"),
             lty = c(1,1,1,1),xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=1)
legend('topright', c("High-risk","Low-risk") , 
       lty = c(1,1), col=c("#E41A1C","BLACK"), bty='n', cex=0.9,lwd = 2, y.intersp=2)
