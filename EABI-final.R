setwd('D:/work/work table/VOD')
library(survival)
library(cmprsk)
library(pec)
library(dplyr)
library(tidyr)
library(DescTools)
library(rms)
library(maxstat)
library(etm)
library(epiDisplay)
load('myEnvironment.RData')
####################bilimax################################
vodh <-tbl_df(read.csv('HDEABIAN1.CSV'))%>%mutate(nrmmn=nrmm-28/30.5,CR=factor(cr),osmn=osm-28/30.5,Bilimaxlog=log2(Bilimax),
                                                  agvhde=aGVHDgrade>0,agvhde3=aGVHDgrade==3)%>%filter(!is.na(eabia),!is.na(age),!is.na(atg),!is.na(ric2rest1),!is.na(Bilimax))

vodhnv<-filter(vodh,vod==0)
vodh28<-filter(vodh,nrmmn>0)
vodh28nv<-filter(vodh28,vod==0)
vodhmax<-filter(vodh,nrmmn>0,vod==0)
boxplot(Bilimax~atg,data=vodhmax)
wilcox.test(Bilimax~atg,data=vodhmax)

coxhdph<-coxph(Surv(nrmmn,nrme)~ bili1vod2 ,data=vodhmax)

logsticr <- glm( eabia~atg,family=binomial(link='logit'),data=vodh28)
logistic.display(logsticr)
cutpoint<-maxgray(vodhmax$nrmmn24m,vodhmax$cr24m, vodhmax$Bilimax, alpha = NA, minprob = 0.1, maxprob = 0.9, event = 2, cens = 0, compevent = 1)



plot.maxgray(cutpoint)
abline(v=3.5)
length(vodhmax$Bilimax)
mstatb <- maxstat.test(Surv(osmn,ose)~ Bilimax ,data=vodhmax,
                       smethod="LogRank",pmethod="exactGauss",
                       abseps=0.01)
plot(mstatb)

coxhdph<-coxph(Surv(nrmmn,nrme)~ Bilimax ,data=vodh28nv)
cox.zph(coxhdph)
########################Bilimax##################################
vodhnv<-filter(vodh, vod==0)

coxhdph<-coxph(Surv(nrmmn,nrme)~ Bilimaxlog ,data=vodb28)
summary(coxhdph)

coxhdph<-coxph(Surv(osmn,ose)~ Bilimaxlog ,data=vodb28)
summary(coxhdph)

coxhdph<-coxph(Surv(aGVHDm,agvhde3)~ Bilimaxlog ,data=vodb)
summary(coxhdph)



#####################eabia in whole cohort day 0 #######################################
vodhsurv <- npsurv(Surv( nrmm, factor(cr))~bili1vod2,data=vodh)
m <-survplot(vodhsurv,state = 2, n.risk=TRUE, xlim=c(0, 120),ylim=c(0, 1),conf.int=0.1,time.inc=20,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.4,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
             col=c('black','red','blue'),
             lty = c(1,1,1,1),xlab='Time After alloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)

vodhsurv <- npsurv(Surv( osm, ose)~bili1vod2,data=vodh)
m <-survplot(vodhsurv, n.risk=TRUE, xlim=c(0, 24),ylim=c(0, 1),conf.int=0.1,time.inc=4,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.4,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
             col=c('black','red','blue'),
             lty = c(1,1,1,1),xlab='Time After alloSCT (Months)', ylab='Survival Probability',cex.xlab=0.9)

#################################day 28#################
vodhsurv <- npsurv(Surv( nrmmn, factor(cr))~eabia,data=vodh28)
m <-survplot(vodhsurv,state = 2, n.risk=TRUE, xlim=c(0, 120),ylim=c(0, 1),conf.int=0.1,time.inc=20,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.4,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
             col=c('black','red','blue'),
             lty = c(1,1,1,1),xlab='Time After Day 28 (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)

vodhsurv <- npsurv(Surv( osmn, ose)~eabia,data=vodh28nv)
m <-survplot(vodhsurv, n.risk=TRUE, xlim=c(0, 24),ylim=c(0, 1),conf.int=0.1,time.inc=4,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.4,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
             col=c('black','red','blue'),
             lty = c(1,1,1,1),xlab='Time After Day 28 (Months)', ylab='Survival Probability',cex.xlab=0.9)

##################p value uni########################################
coxhdph<-coxph(Surv(nrmmn24m,nrme)~ eabia ,data=vodh28)

CIF<-cuminc(ftime =vodh28$nrmmn24m, fstatus = vodh28$cr24m, group =vodh28$eabia)

osfit<- survdiff(Surv( osmn24m, ose24m)~eabia,data=vodh28)

CIFnv<-cuminc(ftime =vodh28nv$nrmmn24m, fstatus = vodh28nv$cr24m, group =vodh28nv$eabia)

osfit<- survdiff(Surv( osmn24m, ose24m)~eabia,data=vodh28nv)

coxos<-coxph(Surv(nrmmn,nrme)~eabia+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data = vodh28nv)
cox.zph(coxos)

coxos<-coxph(Surv(osmn, ose)~eabia+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data = vodh28nv)
summary(coxos)
######################NRM#######################################


hdcoxos<-maxstat_Cox(data=vodh28, maxgray(vodh28$nrmmn,vodh28$CR, vodh28$Bilimax, alpha = NA, minprob = 0.1, maxprob = 0.9, event = 2, cens = 0, compevent = 1),
                     formula=Surv(nrmmn,nrme)~age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1, seed=1119)
summary(hdcoxos$fit)

hdcoxosnv<-maxstat_Cox(data=vodh28nv, maxgray(vodh28$nrmmn,vodh28$CR, vodh28$Bilimax, alpha = NA, minprob = 0.1, maxprob = 0.9, event = 2, cens = 0, compevent = 1),
                     formula=Surv(nrmmn24m,nrme24m)~age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1, seed=1119)
summary(hdcoxosnv$fit)
################OS################
hdcoxosos<-maxstat_Cox(data=vodh28, maxgray(vodh28$nrmmn,vodh28$CR, vodh28$Bilimax, alpha = NA, minprob = 0.1, maxprob = 0.9, event = 2, cens = 0, compevent = 1),
                     formula=Surv(osmn, ose)~age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1, seed=1119)
coxph(Surv(osmn, ose)~eabia+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data=vodh28)
summary(hdcoxosos$fit)

hdcoxososnv<-maxstat_Cox(data=vodh28nv, maxgray(vodh28nv$nrmmn,vodh28nv$CR, vodh28nv$Bilimax, alpha = NA, minprob = 0.1, maxprob = 0.9, event = 2, cens = 0, compevent = 1),
                       formula=Surv(osmn, ose)~age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1, seed=1119)
summary(hdcoxososnv$fit)

hdosraw<-coxph(Surv(osmn, ose)~eabia+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data=vodh28nv)
summary(hdosraw)
#################RELAPSE########################
hdcoxos<-maxstat_Cox(data=vodh28, maxgray(vodh28$nrmmn,vodh28$CR, vodh28$Bilimax, alpha = NA, minprob = 0.1, maxprob = 0.9, event = 2, cens = 0, compevent = 1),
                     formula=Surv(nrmmn,nrme)~age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1, seed=1119)
summary(hdcoxos$fit)
hdcoxos<-maxstat_Cox(data=vodh28, maxgray(vodh28$nrmmn,vodh28$CR, vodh28$Bilimax, alpha = NA, minprob = 0.1, maxprob = 0.9, event = 2, cens = 0, compevent = 1),
                     formula=Surv(nrmmn,nrme)~age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1, seed=1119)
summary(hdcoxos$fit)

hdreraw<-coxph(Surv(nrmmn, relapse)~eabia+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data=vodh28nv)
summary(hdreraw)
##################################TMA+GVHD####################################################

vodhsurv <- npsurv(Surv( nrmmn, factor(cr))~eb1tma2both3,data=vodh28nv)
m <-survplot(vodhsurv,state = 2, n.risk=TRUE, xlim=c(0, 120),ylim=c(0, 1),conf.int=0.1,time.inc=20,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.4,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
             col=c('black','red','blue','brown'),
             lty = c(1,1,1,1),xlab='Time After Day 28 (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)

vodhsurv <- npsurv(Surv( nrmmn, factor(cr))~eb1refrGVHD2both3,data=vodh28nv)
m <-survplot(vodhsurv,state = 2, n.risk=TRUE, xlim=c(0, 120),ylim=c(0, 1),conf.int=0.1,time.inc=20,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.4,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
             col=c('black','red','blue','brown'),
             lty = c(1,1,1,1),xlab='Time After Day 28 (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)

###########################validation cut-off 3.6######################################

vodh28 <-tbl_df(read.csv('HDEABIAN1.CSV'))%>%mutate(nrmmn=nrmm-28/30.5,CR=factor(cr),osmn=osm-28/30.5,osmn24m=ifelse(osmn>24,24,osmn),ose24m=ifelse(osmn>24,0,ose), 
                                                  nrmmn24m=ifelse(nrmmn>24,24,nrmmn),nrme24m=ifelse(nrmm>24,0,nrme),relapse24m=ifelse(nrmmn>24,0,relapse),cr24m=ifelse(nrmmn>24,0,cr)
)%>%filter(!is.na(eabia),!is.na(age),!is.na(atg),!is.na(ric2rest1),!is.na(Bilimax),nrmmn>0)

vodh28nv<-filter(vodh28,vod==0)

hdcoxos<-maxstat_Cox(data=vodh28, maxgray(vodh28$nrmmn,vodh28$CR, vodh28$Bilimax, alpha = NA, minprob = 0.1, maxprob = 0.9, event = 2, cens = 0, compevent = 1),
                     formula=Surv(nrmmn,nrme)~age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1, seed=1119)
summary(hdcoxos$fit)
coxNRM1000bootstrap<-hdcoxos
ls()
save.image(file='myEnvironment.RData')
#####################################berlin#################################

vodb <-tbl_df(read.csv('berlinEABIA.CSV'))%>%mutate(BILI1=(bili1vod2==1)*1,pfs=(cr>0)*1,nrmmn=ifelse(nrmm>=28/30.5,nrmm-28/30.5,0),osmn=ifelse(osm>=28/30.5,osm-28/30.5,0),
                                                    Bilimaxlog=log2(Bilimax),agvhde=aGVHDgrade>0,agvhde3=aGVHDgrade==2)%>%filter(!is.na(age),!is.na(atg),!is.na(ric2rest1),!is.na(BILI1))
vodbnv<-filter(vodb,vod==0)
vodb28<-filter(vodb,nrmmn>0)%>%mutate(off1=0.7894*eabia)
vodb28nv<- filter(vodb28,vod==0)
coxhdph<-coxph(Surv(nrmmn,nrme)~ eabia ,data=vodb)
coxos<-coxph(Surv(nrmmn,nrme)~eabia+offset(off1)+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data = vodb28)
################day 0####################################
vodhsurv <- npsurv(Surv( nrmm, factor(cr))~bili1vod2,data=vodb)
m <-survplot(vodhsurv,state = 2, n.risk=TRUE, xlim=c(0, 50),ylim=c(0, 1),conf.int=0.1,time.inc=10,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.4,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
             col=c('black','red','blue'),
             lty = c(1,1,1,1),xlab='Time After alloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)

CIF<-cuminc(ftime =vodb28nv$nrmmn, fstatus = vodb28nv$cr, group =vodb28nv$eabia)

vodhsurv <- npsurv(Surv( osm, ose)~bili1vod2,data=vodb)
m <-survplot(vodhsurv, n.risk=TRUE, xlim=c(0, 50),ylim=c(0, 1),conf.int=0.1,time.inc=10,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.4,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
             col=c('black','red','blue'),
             lty = c(1,1,1,1),xlab='Time After alloSCT (Months)', ylab='Survival Probability',cex.xlab=0.9)

osfit<- survdiff(Surv( osmn, ose)~eabia,data=vodb28nv)
################day 28####################################
vodhsurv <- npsurv(Surv( nrmmn, factor(cr))~bili1vod2,data=vodb28nv)
m <-survplot(vodhsurv,state = 2, n.risk=TRUE, xlim=c(0, 50),ylim=c(0, 1),conf.int=0.1,time.inc=10,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.4,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
             col=c('black','red','blue'),
             lty = c(1,1,1,1),xlab='Time After Day 28 (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)
vodhsurv <- npsurv(Surv( osmn, ose)~bili1vod2,data=vodb28nv)
m <-survplot(vodhsurv, n.risk=TRUE, xlim=c(0, 50),ylim=c(0, 1),conf.int=0.1,time.inc=10,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.4,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
             col=c('black','red','blue'),
             lty = c(1,1,1,1),xlab='Time After alloSCT (Months)', ylab='Survival Probability',cex.xlab=0.9)

###############################Cox###################################

coxbos<-coxph(Surv(osmn, ose)~eabia+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data = vodb28nv)
summary(coxbos)

coxbnrm<-coxph(Surv(nrmmn, nrme)~eabia,data = vodb28nv)
summary(coxbnrm)

coxbttr<-coxph(Surv(nrmmn, relapse)~eabia+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data = vodb28)
summary(coxbttr)
###############################TIME-DEPENDENT####################################
vodbt<-select(vodb,Nr,eabia,age,sex,mismatch1rest0,aml1rest0,atg,ric2rest1,osm, ose,nrmm, nrme,relapse)

newvodb <- tmerge(data1=vodb, data2=vodb, id=Nr , tstop=nrmm)
newvodb <- tmerge(data1=newvodb, data2=vodb, id=Nr, eabiat=event(eabiam))
newvodba<-mutate(newvodb,eabiat1=1-eabiat*1,nrme=ifelse(nrme==1&tstart==0&eabiat1==1, 0, nrme),
                 ose=ifelse(ose==1&tstart==0&eabiat1==1, 0, ose))

newvodba<-mutate(newvodb,eabiat1=1-eabiat*1,nrme=ifelse(nrme==1&tstart==0&eabiat1==1, 0, nrme),
                 ose=ifelse(ose==1&tstart==0&eabiat1==1, 0, ose),relapse=ifelse(relapse==1&tstart==0&eabiat1==1, 0, relapse))%>%filter(vod==0)

tdcoxnrmUNI<-coxph(Surv(nrmm, nrme)~eabia,data = newvodba)
summary(tdcoxnrmUNI)

tdcoxnrm<-coxph(Surv(nrmm, nrme)~eabia+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data = newvodba)
summary(tdcoxnrm)

tdcoxos<-coxph(Surv(osm, ose)~eabia+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data = newvodba)
summary(tdcoxos)

tdcoxos<-coxph(Surv(osm, ose)~eabia,data = newvodba)
summary(tdcoxos)

tdcoxre<-coxph(Surv(nrmm, relapse)~eabia+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data = newvodba)
summary(tdcoxre)

tdcoxre<-coxph(Surv(nrmm, relapse)~eabia,data = newvodba)
summary(tdcoxre)

#############validation##############################################
vodbva<- mutate(vodb28,offset1=eabia*0.789429)
vodbnvva<- mutate(vodb28nv,off2=eabia*0.74531 )

vacoxnrm<-coxph(Surv(nrmmn, nrme)~eabia+offset(offset1)+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data = vodbva)
summary(vacoxnrm)

vacoxnrm1<-coxph(Surv(nrmmn, nrme)~eabia+offset(off2)+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data = vodbnvva)
summary(vacoxnrm1)
###############################PE validation#######################################

Modelsmultinrmv <- list ("Berlin"=coxph(Surv(nrmmn, nrme)~eabia+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data=vodbva,x=TRUE),
                         "HD offset"=coxph(Surv(nrmmn, nrme)~offset(offset1)+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data=vodbva,x=TRUE))

pehdmultinrmv1<- pec(object=Modelsmultinrmv,
                     formula=Surv(nrmmn, nrme)~age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data=vodbva,reference =FALSE,
                     exact=TRUE,
                     cens.model="cox",
                     splitMethod="Boot632plus",
                     B=1000,
                     verbose=TRUE)
plot(pehdmultinrmv1,col=c('black','red'),lty = c(1,2),legend = FALSE)
legend('topleft', c('Berlin', "HD offset") , 
       lty = c(1,2), col=c('black','red'), bty='n', cex=1,lwd = 2, y.intersp=2)
###########################################################################################
Modelsmultinrmv <- list ("Berlin"=coxph(Surv(nrmmn, nrme)~eabia+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data=vodbnvva,x=TRUE),
                         "HD offset"=coxph(Surv(nrmmn, nrme)~offset(off2)+age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data=vodbnvva,x=TRUE))

pehdmultinrmv1<- pec(object=Modelsmultinrmv,
                     formula=Surv(nrmmn, nrme)~age+sex+mismatch1rest0+aml1rest0+atg+ric2rest1,data=vodbnvva,reference =FALSE,
                     exact=TRUE,
                     cens.model="cox",
                     splitMethod="Boot632plus",
                     B=1000,
                     verbose=TRUE)
plot(pehdmultinrmv1,col=c('black','red'),lty = c(1,2),legend = FALSE)
legend('topleft', c('Berlin', "HD offset") , 
       lty = c(1,2), col=c('black','red'), bty='n', cex=1,lwd = 2, y.intersp=2)

####################logistic ##############################
x<- filter(vodh,!is.na(EASIXpre),statine==1)
logsticr <- glm( eabia~EASIXpre,family='binomial',data=x)
pred.prob <- predict(logsticr,type='response')
brierScore <- mean((pred.prob-x$eabia)^2)
brierScore 
fit <- glm(am~hp+wt,data=mtcars,family='binomial')
pred.prob <- predict(fit,type='response')
brierScore <- mean((pred.prob-mtcars$am)^2)


Brier(list(lm1,lm2),verbose=0,crRatio=1)

logistic.display(logsticr)
wilcox.test(EASIXd0~atg,data=vodbnv)
boxplot(EASIXpre~atg,data=vodbnv)
#####################ROC#######################
data<-vodh
HDbilipre <-roc(vodhnv$eabia, vodhnv$Bilipre,na.rm=TRUE)
summary(HDbilipre)
HDbilipreAUC<-auc(vodhnv$eabia, vodhnv$Bilipre,na.rm=TRUE)
summary(HDbilipreAUC)
HDCHEpre <-ci.auc(vodhnv$eabia, vodhnv$CHEpre,na.rm=TRUE )

a<-rep(0:1,20)
b<-c(rep(0,20),rep(1,20))
sd<-auc(a, b,na.rm=TRUE)
x<-vodh
HDbilipre <-roc(x$eabia, x$Bilipre,na.rm=TRUE)
HDbilipreAUC<-auc(x$eabia, x$Bilipre,na.rm=TRUE)
HDbilipreAUC
HDCHEpre <-ci.auc(x$eabia, x$Bilipre,na.rm=TRUE )
HDCHEpre
roc.test(HDbilipre,sd)

aucall<- function(response,predictor){
  rocoj <-roc(response, predictor,na.rm=TRUE)
  AUC<-auc(response, predictor,na.rm=TRUE)
  ci <-ci.auc(response, predictor,na.rm=TRUE)
out<-list(AUC,ci)
out
}

aucall(vodhnv$eabia,vodhnv$Bilipre)


HDbilipreAUC<-auc(vodhnv$eabia, vodhnv$Bilipre,na.rm=TRUE)
summary(HDbilipreAUC)
HDCHEpre <-ci.auc(vodhnv$eabia, vodhnv$CHEpre,na.rm=TRUE )
HDgGTpre <-roc(vodhnv$eabia, vodhnv$gGTpre,na.rm=TRUE )
HDGPTpre <-roc(vodhnv$eabia, vodhnv$GPTpre,na.rm=TRUE )
###################cut-off bilipre########################
vodhpre<-filter(vodb,!is.na(Bilipre))
cutpoint<-maxgray(vodhpre$nrmmn,vodhpre$cr, vodhpre$Bilipre, alpha = NA, minprob = 0.1, maxprob = 0.9, event = 2, cens = 0, compevent = 1)

plot.maxgray(cutpoint)

#####################AGN2####################
x<- filter(vodh,!is.na(Angiopoietin2pre),statine==1)
logsticr <- glm( eabia~1,family='binomial',data=x)
logistic.display(logsticr)
pred.prob <- predict(logsticr,type='response')
brierScore <- mean((pred.prob-x$eabia)^2)
brierScore 
length(pred.prob)
######################statine################################
x<- filter(vodh28nv,statine==1)

vodhsurv <- npsurv(Surv( nrmmn, factor(cr))~eabia,data=x)
m <-survplot(vodhsurv,state = 2, n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.4,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
             col=c('black','red','blue'),
             lty = c(1,1,1,1),xlab='Time After Day 28 (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)


vodhsurv <- npsurv(Surv( osmn, ose)~eabia,data=x)
m <-survplot(vodhsurv, n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=FALSE,lwd = 2,
             y.n.risk=-0.4,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
             col=c('black','red','blue'),
             lty = c(1,1,1,1),xlab='Time After Day 28 (Months)', ylab='Survival Probability',cex.xlab=0.9)



CIF<-cuminc(ftime =vodh28$nrmmn24m, fstatus = vodh28$cr24m, group =vodh28$eabia)

osfit<- survdiff(Surv( osmn24m, ose24m)~eabia,data=vodh28)

##############HD+berlin######################
vodhc1<- filter(vodh28,eabia==1,statine==0)%>%dplyr::select(nrmmn,cr,vod)%>%mutate(group=1)
vodhc2<- filter(vodh28,eabia==1,statine==1)%>%dplyr::select(nrmmn,cr,vod)%>%mutate(group=2)
vodhc3<- filter(vodh28,eabia==0)%>%dplyr::select(nrmmn,cr,vod)%>%mutate(group=4)
vodbc1<-filter(vodb28,eabia==1)%>%dplyr::select(nrmmn,cr,vod)%>%mutate(group=3)
vodbc2<-filter(vodb28,eabia==0)%>%dplyr::select(nrmmn,cr,vod)%>%mutate(group=4)
vodcomb<- rbind(vodhc1,vodhc2,vodhc3)
vodcombnv<-filter(vodcomb,vod==0)
vodh28s<-mutate(vodh28,es=(eabia+1)*(statine+2))

vodhsurv <- npsurv(Surv( nrmmn, factor(cr))~es,data=vodh28s)
m <-survplot(vodhsurv,state = 2, n.risk=TRUE, xlim=c(0, 24),ylim=c(0, 1),conf.int=0.1,time.inc=4,label.curves=TRUE,lwd = 2,
             y.n.risk=-0.4,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
             col=c('black','blue','red','brown'),
             lty = c(1,1,1,1),xlab='Time After Day 28 (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)




