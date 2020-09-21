setwd('D:/work/work table/ADMA')
library(survival)
library(cmprsk)
library(pec)
library(dplyr)
library(DescTools)
library(rms)
library(maxstat)
library(etm)
adma <- tbl_df(read.csv('HEADMA.CSV'))
adma <- mutate(adma, admalog2=log2(adma$ADMA))

a<-prodlim(Surv(timeTPLuntildeadLFU,Ose)~1,data=adma)
summary(a)
adma <- mutate(adma, pfs=CR>0)
adma$CR <- as.factor(adma$CR)
os.cox <- coxph(Surv(timeTPLuntildeadLFU, Ose)~admalog2+age+RecipientSex+noMismatch1Rest0+dosex
                 +MAC.Aplist1RICist2+GratwohlSCORE+strata(center),data=adma)
os.cox <- npsurv(Surv(timeTPLrelapse.NRM , CR)~1,data=adma)
?npsurv
summary(a, times=c(12))
m <-survplot(os.cox, state = 2, n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='black',sep.n.risk=0.06, cex.axis=2,
         col=c('#969696','#969696','black','black'),
         lty = c(1,2,2,1),xlab='Time After AlloSCT (Months)', ylab='survival probability',cex.xlab=0.9)

k <- cox.zph(os.cox)
k
plot(k[1],lwd=2,xlab='Time after alloSCT (months)')
abline(0,0, col=2,lwd=2)
abline(h= os.cox$coef[1], col=3, lwd=2, lty=2)
legend('topright', c("admalog2 effect" , 'CI for admalog2effect',"no effect","average effect" ) , 
       lty = c(1,2,1,2), col=c(' black',' black','red','green' ), bty='n', cex=1,lwd = 4, y.intersp=1)
admad <- survSplit(Surv(timeTPLuntildeadLFU, Ose) ~ ., data= adma, cut=c(12,36),
                   episode= "tgroup", id="nr")
hdadmad <- filter(admad, center=='H')
enadmad <- filter(admad, center=='E')
t1 <- filter(admad, tgroup==1)
OSCUT<- coxph(Surv(tstart,timeTPLuntildeadLFU, Ose)~admalog2:strata(tgroup)+age+RecipientSex+noMismatch1Rest0+dosex
                  +GratwohlSCORE+MAC.Aplist1RICist2,data=admad)
mstat <- maxstat.test(Surv(timeTPLuntildeadLFU, Ose) ~ admalog2, data=t1,
                      smethod="LogRank", pmethod="exactGauss",
                      abseps=0.01)
OSCUT<- coxph(Surv(tstart,timeTPLuntildeadLFU, Ose)~admalog2+age+RecipientSex+noMismatch1Rest0+dosex
              +GratwohlSCORE+MAC.Aplist1RICist2,data=t1)
plot(mstat)
summary(OSCUT)
admad1 <-filter(admad, tgroup==1)
Modelsos <- list("Cox.ADMA"=coxph(Surv(tstart,timeTPLuntildeadLFU, Ose)~admalog2+strata(center),data=admad1))
PredError.632plus <- pec(object=Modelsos,
                         formula=Surv(timeTPLuntildeadLFU, Ose)~admalog2,data=admad1, reference = TRUE,
                         exact=TRUE,
                         cens.model="marginal",
                         splitMethod="non",
                         B=100,
                         verbose=TRUE)
plot(PredError.632plus,xlim = c(0,11.9),legend=FALSE,lty = c(1,2,1,2),lwd = 2,xlab = 'Time After AlloSCT')
legend('topleft', c("Reference" ,"ADMA") , 
       lty = c(1,2), col=c(' black','red' ), bty='n', cex=1,lwd = 2, y.intersp=2)

summary(OSCUT)
cox.zph(OSCUT)
pfs.cox <- coxph(Surv(timeTPLrelapse.NRM, pfs)~admalog2+age+RecipientSex+noMismatch1Rest0+dosex
                 +MAC.Aplist1RICist2+GratwohlSCORE+strata(center),data=adma)
k2 <- cox.zph(pfs.cox)
k2
plot(k2[1],lwd=2,xlab='Time after alloSCT (months)')
abline(0,0, col=2,lwd=2)
abline(h= pfs.cox$coef[1], col=3, lwd=2, lty=2)
admapfs <- survSplit(Surv(timeTPLrelapse.NRM, pfs) ~ ., data= adma, cut=c(12,36),
                   episode= "tgroup", id="nr")
PFS1 <- mutate(admapfs,tgroup==1)
hdadmapfs <- filter(admapfs, center=='H')
enadmapfs<- filter(admapfs, center=='E')
PFSCUT<- coxph(Surv(tstart,timeTPLrelapse.NRM,pfs)~admalog2:strata(tgroup)+age+RecipientSex+noMismatch1Rest0+dosex
              +GratwohlSCORE+MAC.Aplist1RICist2,data=admapfs)
summary(PFSCUT)
cox.zph(PFSCUT)
admadpfs1 <-filter(admapfs, tgroup==1) 
Modelspfs <- list("Cox.ADMA"=coxph(Surv(timeTPLuntildeadLFU, Ose)~admalog2+strata(center),data=admadpfs1))
PredError.632plus1 <- pec(object=Modelspfs,
                         formula=Surv(timeTPLuntildeadLFU, Ose)~admalog2,data=admadpfs1, reference = TRUE,
                         exact=TRUE,
                         cens.model="marginal",
                         splitMethod="non",
                         B=100,
                         verbose=TRUE)
plot(PredError.632plus1,xlim = c(0,11.9),legend=FALSE,lty = c(1,2,1,2),lwd = 2,xlab = 'Time After AlloSCT')
legend('topleft', c("Reference" ,"ADMA") , 
       lty = c(1,2), col=c(' black','red' ), bty='n', cex=1,lwd = 2, y.intersp=2)
NRM
admapfs
sum(adma$NRMe)
sum(admapfs$NRMe)
admanrm <- survSplit(Surv(timeTPLrelapse.NRM, NRMe) ~ ., data= adma, cut=c(12,36),
                     episode= "tgroup", id="nr")
admanrm1 <- filter(admanrm, tgroup==1)
hdadmanrm <- filter(admanrm, center=='H')
enadmanrm<- filter(admanrm, center=='E')
NRM<- coxph(Surv(tstart,timeTPLrelapse.NRM,NRMe)~admalog2:strata(tgroup)+age+RecipientSex+noMismatch1Rest0+dosex
               +GratwohlSCORE+MAC.Aplist1RICist2,data=admanrm)

nrm1 <-cuminc(admanrm$timeTPLrelapse.NRM, admanrm$NRMe, rho=0, cencode=0,
        na.action=na.omit)
timepoints(nrm1, 12)
plot(nrm1)
summary(NRM)
cox.zph(NRM)
admarelapse <- survSplit(Surv(timeTPLrelapse.NRM, Relapse) ~ ., data= adma, cut=c(12,36),
                     episode= "tgroup", id="nr")
hdadmarelapse <- filter(admarelapse, center=='H')
enadmarelapse<- filter(admarelapse, center=='E')
relaspe<- coxph(Surv(tstart,timeTPLrelapse.NRM,Relapse)~admalog2:strata(tgroup)+age+RecipientSex+noMismatch1Rest0+dosex
            +GratwohlSCORE+MAC.Aplist1RICist2,data=admarelapse)
summary(relaspe)
cox.zph(relaspe)
##################OSA####################
os.cox1 <- coxph(Surv(Osafter_aGVHD, Ose)~admalog2+age+RecipientSex+noMismatch1Rest0+dosex
                +MAC.Aplist1RICist2+GratwohlSCORE+strata(center),data=adma)
k <- cox.zph(os.cox1)
k
plot(k[1],lwd=2,xlab='Time after aGVHD (months)')
abline(0,0, col=2,lwd=2)
abline(h= os.cox1$coef[1], col=3, lwd=2, lty=2)
admaosa <- filter(adma,!is.na(Osafter_aGVHD))
admaosa <- filter(admaosa, admaosa$Osafter_aGVHD>0)


admaosa<- survSplit(Surv(Osafter_aGVHD, Ose) ~ ., data= admaosa, cut=c(12,36),
                               episode= "tgroup", id="nr")
osa1 <-mutate(admaosa,tgroup==1)
hdadmaosa <- filter(admaosa, center=='H')
enadmaosa<- filter(admaosa, center=='E')
osa<- coxph(Surv(tstart,Osafter_aGVHD, Ose)~admalog2:strata(tgroup)+age+RecipientSex+noMismatch1Rest0+dosex
                +GratwohlSCORE+MAC.Aplist1RICist2,data=admaosa)
summary(osa)
cox.zph(osa)
admaosa1 <-filter(admaosa, tgroup==1) 
Modelsosa <- list("Cox.ADMA"=coxph(Surv(Osafter_aGVHD, Ose)~admalog2+strata(center),data=admaosa1))
PredError.632plus1 <- pec(object=Modelsosa,
                          formula=Surv(Osafter_aGVHD, Ose)~admalog2,data=admaosa1, reference = TRUE,
                          exact=TRUE,
                          cens.model="marginal",
                          splitMethod="non",
                          B=100,
                          verbose=TRUE)
plot(PredError.632plus1,xlim = c(0,11.9),legend=FALSE,lty = c(1,2,1,2),lwd = 2,xlab = 'Time After aGVHD')
legend('topleft', c("Reference" ,"ADMA") , 
       lty = c(1,2), col=c(' black','red' ), bty='n', cex=1,lwd = 2, y.intersp=2)
#################NRMA#####################
pfsa.cox1 <- coxph(Surv(NRMafter_aGVHD, pfsa)~admalog2+age+RecipientSex+noMismatch1Rest0+dosex
                 +MAC.Aplist1RICist2+GratwohlSCORE+strata(center),data=adma)
k1 <- cox.zph(pfsa.cox1)
k1
plot(k1[1],lwd=2,xlab='Time after aGVHD (months)')
abline(0,0, col=2,lwd=2)
abline(h= pfsa.cox1$coef[1], col=3, lwd=2, lty=2)
adma1 <- filter(adma,!is.na(NRMafter_aGVHD))
adma2 <- mutate(adma1, nrma=cra==2)
adma3 <- filter(adma2, adma2$NRMafter_aGVHD>0)
adma$nrma
admanrma<- survSplit(Surv(NRMafter_aGVHD, nrma) ~ ., data= adma3, cut=c(12,36),
                    episode= "tgroup", id="nr")
hdadmanrma <- filter(admanrma, center=='H')
enadmanrma<- filter(admanrma, center=='E')
nrma<- coxph(Surv(tstart,NRMafter_aGVHD, nrma)~admalog2:strata(tgroup)+age+RecipientSex+noMismatch1Rest0+dosex
            +GratwohlSCORE+MAC.Aplist1RICist2,data=admanrma)
summary(nrma)
cox.zph(nrma)

adma1 <- filter(adma, adma$NRMafter_aGVHD>0)
admapfsa <- mutate(adma1,pfsa=cra>0)


admapfsa<- survSplit(Surv(NRMafter_aGVHD, pfsa) ~ ., data= admapfsa, cut=c(12,36),
                     episode= "tgroup", id="nr")
admapfsa1 <- mutate(admapfsa,tgroup==1)
hdadmapfsa <- filter(admapfsa, center=='H')
enadmapfsa<- filter(admapfsa, center=='E')
pfsa<- coxph(Surv(tstart,NRMafter_aGVHD, pfsa)~admalog2:strata(tgroup)+age+RecipientSex+noMismatch1Rest0+dosex
             +GratwohlSCORE+MAC.Aplist1RICist2,data=admapfsa)
summary(pfsa)
cox.zph(pfsa)
admapfsa1 <-filter(admapfsa, tgroup==1) 
Modelpfsa <- list("Cox.ADMA"=coxph(Surv(NRMafter_aGVHD, pfsa)~admalog2+strata(center),data=admapfsa1))
PredError.632plus1 <- pec(object=Modelpfsa,
                          formula=Surv(NRMafter_aGVHD, pfsa)~admalog2,data=admapfsa1, reference = TRUE,
                          exact=TRUE,
                          cens.model="marginal",
                          splitMethod="non",
                          B=100,
                          verbose=TRUE)
plot(PredError.632plus1,xlim = c(0,11.9),legend=FALSE,lty = c(1,2,1,2),lwd = 2,xlab = 'Time After aGVHD')
legend('topleft', c("Reference" ,"ADMA") , 
       lty = c(1,2), col=c(' black','red' ), bty='n', cex=1,lwd = 2, y.intersp=2)

pfsa<- coxph(Surv(tstart,NRMafter_aGVHD, pfsa)~admalog2:strata(tgroup)+age+RecipientSex+noMismatch1Rest0+dosex
             +GratwohlSCORE+MAC.Aplist1RICist2,data=admapfsa)
#####################SURV PLOTS###################################
adma <- tbl_df(read.csv('HEADMA.CSV'))
adma <- mutate(adma, admalog2=log2(adma$ADMA))
adma <- mutate(adma, pfs=CR>0)
adma$CR <- as.factor(adma$CR)
adma$cra <- as.factor(adma$cra)
adma$Q <- CutQ(adma$admalog2 , breaks = quantile(adma$admalog2, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)

summary(adma$admalog2)
hist(adma$admalog2,c)
admad1
fitos <- npsurv(Surv(timeTPLuntildeadLFU, Ose)~Q,data=adma)
par(mgp=c(1.5, 1, 0))
survplot(fitos,  n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=TRUE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('#969696','#969696','black','black'),
         lty = c(1,2,2,1),xlab='Time After AlloSCT (Months)', ylab='survival probability',cex.xlab=0.9)
legend('bottomleft', c("Q1","Q2",'Q3', "Q4") , 
       lty = c(1,2,2,1), col=c('#969696','#969696','black','black'), bty='n', cex=0.9,lwd = 2, y.intersp=2)

admapfs <- survSplit(Surv(timeTPLrelapse.NRM, pfs) ~ ., data= adma, cut=c(12,36),
                     episode= "tgroup", id="nr")
admadpfs1 <-filter(admapfs, tgroup==1) 
admadpfs1 <- mutate(admadpfs1, ros=Ose*Relapse)
fitos <- npsurv(Surv(timeTPLrelapse.NRM, CR)~Q,data=adma)
par(mgp=c(1.5, 1, 0))
survplot(fitos,state=2,  n.risk=TRUE, xlim=c(0, 100),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=TRUE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('#969696','#969696','black','black'),
         lty = c(1,2,2,1),xlab='Time After AlloSCT (Months)', ylab='survival probability',cex.xlab=0.9)
legend('bottomleft', c("Q1","Q2",'Q3', "Q4") , 
       lty = c(1,2,2,1), col=c('#969696','#969696','black','black'), bty='n', cex=0.9,lwd = 2, y.intersp=2)
admanrm <- survSplit(Surv(timeTPLrelapse.NRM, CR) ~ ., data= adma, cut=c(12.1,36),
                     episode= "tgroup", id="nr")
admanrm1 <-filter(admanrm, tgroup==1) 

summary(fitCI,times = 12)
admanrm1$CR
par(mgp=c(2, 1, 0))
summary(fitos)
survplot(fitos,  state = '1',n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='black',sep.n.risk=0.06, cex.axis=2,
         col=c('#969696','#969696','black','black'),
         lty = c(1,2,2,1),xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)
legend('topleft', c("Q1","Q2",'Q3', "Q4") , 
       lty = c(1,2,2,1), col=c('#969696','#969696','black','black'), bty='n', cex=0.9,lwd = 2, y.intersp=2)

admaosa <- filter(adma,!is.na(Osafter_aGVHD))
admaosa <- filter(admaosa, admaosa$Osafter_aGVHD>0)
admaosa<- survSplit(Surv(Osafter_aGVHD, Ose) ~ ., data= admaosa, cut=c(12,36),
                    episode= "tgroup", id="nr")
admaosa1 <-filter(admaosa, tgroup==1) 
fitos <- npsurv(Surv(Osafter_aGVHD, Ose)~Q,data=admaosa1)
par(mgp=c(1.5, 1, 0))
survplot(fitos,  n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('#969696','#969696','black','black'),
         lty = c(1,2,2,1),xlab='Time After aGVHD (Months)', ylab='survival probability',cex.xlab=0.9)
legend('bottomleft', c("Q1","Q2",'Q3', "Q4") , 
       lty = c(1,2,2,1), col=c('#969696','#969696','black','black'), bty='n', cex=0.9,lwd = 2, y.intersp=2)

adma <- filter(adma, adma$NRMafter_aGVHD>0)
admapfsa <- mutate(adma,pfsa=cra>0)

admapfsa<- survSplit(Surv(NRMafter_aGVHD, pfsa) ~ ., data= admapfsa, cut=c(12,36),
                     episode= "tgroup", id="nr")
admapfsa1 <-filter(admapfsa, tgroup==1) 
fitos <- npsurv(Surv(NRMafter_aGVHD, pfsa)~Q,data=admapfsa1)
survplot(fitos,  n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('#969696','#969696','black','black'),
         lty = c(1,2,2,1),xlab='Time After aGVHD (Months)', ylab='survival probability',cex.xlab=0.9)
legend('bottomleft', c("Q1","Q2",'Q3', "Q4") , 
       lty = c(1,2,2,1), col=c('#969696','#969696','black','black'), bty='n', cex=0.9,lwd = 2, y.intersp=2)
adma <- filter(adma,!is.na(NRMafter_aGVHD))

adma <- filter(adma, adma$NRMafter_aGVHD>0)
admanrma<- survSplit(Surv(NRMafter_aGVHD, cra) ~ ., data= adma, cut=c(12,36),
                     episode= "tgroup", id="nr")
admanrma1 <-filter(admanrma, tgroup==1) 
fitos <- npsurv(Surv(NRMafter_aGVHD, cra)~Q,data=admanrma1)
survplot(fitos,  state = '2',n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('#969696','#969696','black','black'),
         lty = c(1,2,2,1),xlab='Time After aGVHD (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)
legend('topleft', c("Q1","Q2",'Q3', "Q4") , 
       lty = c(1,2,2,1), col=c('#969696','#969696','black','black'), bty='n', cex=0.9,lwd = 2, y.intersp=2)

##############iNOS##################
#########rs1137933
admad1<-  mutate(admad1,newa=Q=='Q4',rs1137933=rs1137933_iNOS_SNP)
admad1
admad933<- filter(admad1, !is.na(admad1$rs1137933_iNOS_SNP))

admad933h <- filter(admad933,newa==TRUE)
admad933l <- filter(admad933,newa==FALSE)
admad933h
fitos <- npsurv(Surv(timeTPLuntildeadLFU, Ose)~rs1137933,admad933h)
fitos <- npsurv(Surv(timeTPLuntildeadLFU, Ose)~rs1137933,admad933l)
par(mgp=c(1.5, 1, 0))
survplot(fitos,  n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('black','red','blue'),lty=1,
         xlab='Time After AlloSCT (Months)', ylab='survival probability',cex.xlab=0.9)

legend('bottomleft', c("GG","AG",'AA') , 
       lty = 1, col=c('black','red','blue'), bty='n', cex=0.9,lwd = 2, y.intersp=2)

adma$CR<- as.factor(adma$CR)
admanrm <- survSplit(Surv(timeTPLrelapse.NRM, CR) ~ ., data= adma, cut=c(12.1,36),
                     episode= "tgroup", id="nr")
admanrm1 <- filter(admanrm,tgroup==1)
admanrm1<-  mutate(admanrm1,newa=Q=='Q4',rs1137933=rs1137933_iNOS_SNP)
admad933<- filter(admanrm1, !is.na(admanrm1$rs1137933_iNOS_SNP))
admad933$rs1137933
admad933h <- filter(admad933,newa==TRUE)
admad933h
admad933l <- filter(admad933,newa==FALSE)
admad933h

fitos <- npsurv(Surv(timeTPLrelapse.NRM, CR)~rs1137933,admad933h)
fitos <- npsurv(Surv(timeTPLrelapse.NRM, CR)~rs1137933,admad933l)
par(mgp=c(1.5, 1, 0))
survplot(fitos, state = '2' ,n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('black','red','blue'),lty=1,
         xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)

legend('topleft', c("GG","AG",'AA'), 
       lty = 1, col=c('black','red','blue'), bty='n', cex=0.9,lwd = 2, y.intersp=2)



admad

admad1<-  mutate(admad1,newa=Q=='Q4',rs2297518=4-rs2297518_iNOS_SNP)
admad1
admad518<- filter(admad1, !is.na(admad1$rs2297518_iNOS_SNP))

admad518h <- filter(admad518,newa==TRUE)
admad518l <- filter(admad518,newa==FALSE)
admad933h
fitos <- npsurv(Surv(timeTPLuntildeadLFU, Ose)~rs2297518,admad518h)
fitos <- npsurv(Surv(timeTPLuntildeadLFU, Ose)~rs2297518,admad518l)
par(mgp=c(1.5, 1, 0))
survplot(fitos,  n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('black','red','blue'),lty=1,
         xlab='Time After AlloSCT (Months)', ylab='survival probability',cex.xlab=0.9)

legend('bottomleft', c("GG","AG",'AA') , 
       lty = 1, col=c('black','red','blue'), bty='n', cex=0.9,lwd = 2, y.intersp=2)



adma$CR<- as.factor(adma$CR)
admanrm <- survSplit(Surv(timeTPLrelapse.NRM, CR) ~ ., data= adma, cut=c(12.1,36),
                     episode= "tgroup", id="nr")
admanrm1 <- filter(admanrm,tgroup==1)
admanrm1<-  mutate(admanrm1,newa=Q=='Q4',rs2297518=4-rs2297518_iNOS_SNP)
admad518<- filter(admanrm1, !is.na(admanrm1$rs2297518_iNOS_SNP))
admad518$rs1137933
admad518h <- filter(admad518,newa==TRUE)
admad518h
admad518l <- filter(admad518,newa==FALSE)
admad518h

fitos <- npsurv(Surv(timeTPLrelapse.NRM, CR)~rs2297518,admad518h)
fitos <- npsurv(Surv(timeTPLrelapse.NRM, CR)~rs2297518,admad518l)
par(mgp=c(1.5, 1, 0))
survplot(fitos, state = '2' ,n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('black','red','blue'),lty=1,
         xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)

legend('topleft', c("GG","AG",'AA'), 
       lty = 1, col=c('black','red','blue'), bty='n', cex=0.9,lwd = 2, y.intersp=2)


admad1<-  mutate(admad1,newa=Q=='Q4',rs2779248=4-rs2779248_iNOS_SNP)
admad1
admad248<- filter(admad1, !is.na(admad1$rs2779248_iNOS_SNP))

admad248h <- filter(admad248,newa==TRUE)
admad248l <- filter(admad248,newa==FALSE)

fitos <- npsurv(Surv(timeTPLuntildeadLFU, Ose)~rs2779248,admad248h)
fitos <- npsurv(Surv(timeTPLuntildeadLFU, Ose)~rs2779248,admad248l)
par(mgp=c(1.5, 1, 0))
survplot(fitos,  n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('black','red','blue'),lty=1,
         xlab='Time After AlloSCT (Months)', ylab='survival probability',cex.xlab=0.9)

legend('bottomleft', c("TT","CT",'CC') , 
       lty = 1, col=c('black','red','blue'), bty='n', cex=0.9,lwd = 2, y.intersp=2)



adma$CR<- as.factor(adma$CR)
admanrm <- survSplit(Surv(timeTPLrelapse.NRM, CR) ~ ., data= adma, cut=c(12.1,36),
                     episode= "tgroup", id="nr")
admanrm1 <- filter(admanrm,tgroup==1)
admanrm1<-  mutate(admanrm1,newa=Q=='Q4',rs2779248=rs2779248_iNOS_SNP)
admad248<- filter(admanrm1, !is.na(admanrm1$rs2779248_iNOS_SNP))
admad248$rs1137933
admad248h <- filter(admad248,newa==TRUE)
admad248h
admad248l <- filter(admad248,newa==FALSE)
admad248h

fitos <- npsurv(Surv(timeTPLrelapse.NRM, CR)~rs2779248,admad248h)
fitos <- npsurv(Surv(timeTPLrelapse.NRM, CR)~rs2779248,admad248l)
par(mgp=c(1.5, 1, 0))
survplot(fitos, state = '2' ,n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('black','red','blue'),lty=1,
         xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)

legend('topleft', c("TT","CT",'CC'), 
       lty = 1, col=c('black','red','blue'), bty='n', cex=0.9,lwd = 2, y.intersp=2)
###########likelihood ratio test#######
OSCUT1<- coxph(Surv(tstart,timeTPLuntildeadLFU, Ose)~admalog2+age+RecipientSex+noMismatch1Rest0+dosex
              +GratwohlSCORE+MAC.Aplist1RICist2,data=t1)
OSCUT2<- coxph(Surv(tstart,timeTPLuntildeadLFU, Ose)~age+RecipientSex+noMismatch1Rest0+dosex
              +GratwohlSCORE+MAC.Aplist1RICist2,data=t1)
anova(OSCUT2,OSCUT1)

PFSCUT1<- coxph(Surv(tstart,timeTPLrelapse.NRM,pfs)~admalog2+age+RecipientSex+noMismatch1Rest0+dosex
               +GratwohlSCORE+MAC.Aplist1RICist2,data=PFS1)
PFSCUT2<- coxph(Surv(tstart,timeTPLrelapse.NRM,pfs)~age+RecipientSex+noMismatch1Rest0+dosex
                +GratwohlSCORE+MAC.Aplist1RICist2,data=PFS1)
anova(PFSCUT1,PFSCUT2)

osaa1<- coxph(Surv(tstart,Osafter_aGVHD, Ose)~admalog2+age+RecipientSex+noMismatch1Rest0+dosex
            +GratwohlSCORE+MAC.Aplist1RICist2,data=osa1)
osaa2<- coxph(Surv(tstart,Osafter_aGVHD, Ose)~age+RecipientSex+noMismatch1Rest0+dosex
             +GratwohlSCORE+MAC.Aplist1RICist2,data=osa1)
anova(osaa1,osaa2)

pfsa1<- coxph(Surv(tstart,NRMafter_aGVHD, pfsa)~admalog2+age+RecipientSex+noMismatch1Rest0+dosex
             +GratwohlSCORE+MAC.Aplist1RICist2,data=admapfsa)
pfsa2<- coxph(Surv(tstart,NRMafter_aGVHD, pfsa)~age+RecipientSex+noMismatch1Rest0+dosex
              +GratwohlSCORE+MAC.Aplist1RICist2,data=admapfsa)
anova(pfsa1,pfsa2)
