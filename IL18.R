setwd('D:/work/work table/IL18')
library(survival)
library(cmprsk)
library(pec)
library(dplyr)
library(DescTools)
library(rms)
library(maxstat)
library(etm)
library(pROC)
library(ggsci)
############BOXPLOT####################
ilhdbox <- tbl_df(read.csv('il18boxplot.CSV'))
ilhdbox$centre<- ordered(ilhdbox$centre, levels=c('N','H','E'))
ilhdboxn<-filter(ilhdbox,centre=='E')
summary(ilhdboxn$freeIL18pre)
ggplot(data =filter(ilhdbox,!is.na(freeIL18pre)), aes(x=factor(centre), y=freeIL18pre)) + geom_boxplot()+theme_bw()  
ggplot(data =filter(ilhdbox,!is.na(IL18BPapre)), aes(x=factor(centre), y=IL18BPapre)) + geom_boxplot()+theme_bw()  
ggplot(data =filter(ilhdbox,!is.na(IL18pre)), aes(x=factor(centre), y=IL18pre)) + geom_boxplot()+theme_bw() 


ilhd <- mutate(tbl_df(read.csv('IL18FINAL.CSV')),il18prelog=log2(freeIL18pre),age=Age,cr1=as.factor(cr),freeIL18d28=as.numeric(freeIL18d28),
               Gratwohl2vsrest=Gratwohl2vsrest==2, mismatch3rest1=mismatch3rest1==3, donorsex=donorsex==2,sex=sex==2,myeloid1lymphoid2=myeloid1lymphoid2==2,
               atg=atg==1,ric2rest1=ric2rest1==2)

summary(ilhd$il18prelog)
ilhd$con13
ilhdnn<- filter(ilhd, !is.na(freeIL18pre), !is.na(freeIL18d28))
con13n==1|con13n==5
summary(ilhdnn$freeIL18pre)
summary(ilhdnn$freeIL18d28)
ilhds<- select(ilhd,freeIL18pre, freeIL18d28)%>%gather(timep, il18)
ilhds$timep<-ordered(ilhds$timep,levels=c('freeIL18pre','freeIL18d28'))
boxplot(il18~timep,data=ilhds)
summary(factor(ilhd$con13n))
ilhdwt<- filter(ilhd,!is.na(freeIL18pre),con13n==1)
wilcox.test(freeIL18pre~con13n,data=ilhdwt)
summary(factor(ilhdwt$con13n))
summary(ilhdwt$freeIL18pre)

boxplot(freeIL18pre~con13n,data=ilhd)
wilcox.test(freeIL18pre~con13,data=ilhdwt) 
wilcox.test(freeIL18d28~con13,data=ilhdwt) 
summary(freeIL18pre~con13,data=ilhd)

ilhd$Qpre <- CutQ(ilhd $freeIL18pre, breaks = quantile(ilhd$freeIL18pre, seq(1, 0, by = -0.25), na.rm = TRUE), labels = NULL, na.rm = FALSE)

ilhd$Qpre <- CutQ(ilhd $freeIL18d28, breaks = quantile(ilhd$freeIL18d28, seq(1, 0, by = -0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)


head(ilhd)
boxplot(freeIL18pre~con13,data=ilhd) 

summary(a)
a<-ggplot(data =filter(ilhd,!is.na(freeIL18pre)), aes(x=factor(con13n), y=freeIL18pre)) + geom_boxplot()+theme_bw()

ggplot(data =filter(ilhd,!is.na(freeIL18d28)), aes(x=factor(con1), y=freeIL18d28)) + geom_boxplot()+theme_bw()
ggplot(data =filter(ilhd,!is.na(freeIL18d28)), aes(x=factor(con13), y=freeIL18d28)) + geom_boxplot()+theme_bw()

########Essen####################
ggplot(data =filter(ilen,!is.na(freeIL18pre)), aes(x=factor(con1y13), y=freeIL18pre)) + geom_boxplot()+theme_bw()
summary(ilen$freeIL18pre~ilen$con1y13)

ilenwt<- filter(ilen,!is.na(freeIL18pre),con1y13==1)
boxplot(freeIL18pre~con1y13,data=ilenwt)
wilcox.test(freeIL18pre~con1y13,data=ilenwt) 
summary(ilenwt$freeIL18pre)
################cause of death bar plot#########
library(reshape)
library(scales)
library(ggsci)
datm <- melt(cbind(dat, ind = rownames(dat)), id.vars = c('ind'))
ilhd1<-filter(ilhd,!is.na(con13n),nrme1y==1,!is.na(freeIL18pre))
summary(ilhd1$freeIL18pre)
summary(factor(ilhd1$con13n))
ilhd1$Qfreeil18pre1<-CutQ(ilhd1$il18prelog, breaks=quantile(ilhd1$il18prelog, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhd1$con13n <- ordered(ilhd1$con13n ,levels=c(5,3,4,1))
ilhd1bar<- mutate(ilhd1, q14=Qfreeil18pre1!='Q1',count=1)%>%filter(!is.na(Qfreeil18pre1),Qfreeil18pre1=='Q1'|Qfreeil18pre1=='Q4')

percentilhd1bar <- ilhd1bar %>% group_by(Qfreeil18pre1) %>% count(con13n) %>%mutate(ratio=scales::percent(n/sum(n)))

ggplot(ilhd1bar,aes(x = Qfreeil18pre1, fill = as.factor(con13n))) + geom_bar(position = "fill",width = 0.5)+scale_fill_jco()+
scale_y_continuous(labels = percent_format())+ geom_text(data=percentilhd1bar, aes(y=n,label=ratio),
position=position_fill(vjust=0.5))

ilhd$freeIL18pre
date1<- filter(ilhd,con1==4)
summary(date1$freeIL18pre)
date1<- filter(ilhd,con1==4)
summary(date1$freeIL18d28)
#######BARplot ESSEN##############################
ilen <- mutate(tbl_df(read.csv('IL18EN.CSV')),il18prelog=log2(freeIL18pre),cr1=as.factor(cr),
               Gratwohl2vsrest=Gratwohl2vsrest==2, mismatch3rest1=mismatch3rest1==3, donorsex=donorsex==2,
               sex=sex==2,myeloid1lymphoid2=myeloid1lymphoid2==2,IL18BPaprelog=log2(IL18BPapre), allil18prelog=log2(IL18pre),
               atg=atg==1,ric2rest1=ric2rest1==2)
ilen1<-filter(ilen,!is.na(con13),nrme1y==1)
dateilen1<- filter(ilen0,con1y==0)
summary(ilen1$freeIL18pre)
ilen1$Qfreeil18pre1<-CutQ(ilen1$il18prelog, breaks=quantile(ilen1$il18prelog, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilen1$con13 <- ordered(ilen1$con13 ,levels=c(5,3,4,1))
ilen1bar<- mutate(ilen1, q14=Qfreeil18pre1!='Q1',count=1)%>%filter(!is.na(q14),Qfreeil18pre1=='Q1'|Qfreeil18pre1=='Q4')

ilen1bar<- mutate(ilen1, q14=freeIL18pre>=903.60,count=1)%>%filter(!is.na(q14))

percentilhd1bar <- ilen1bar %>% group_by(q14) %>% count(con13) %>%
  mutate(ratio=scales::percent(n/sum(n)))

ggplot(ilen1bar,aes(x = q14, fill = as.factor(con13))) + 
  geom_bar(position = "fill",width = 0.5)+scale_fill_jco()+
  scale_y_continuous(labels = percent_format())+ geom_text(data=percentilhd1bar, aes(y=n,label=ratio),
                                                           position=position_fill(vjust=0.5))
###################median follow up and incidences#####################

fit <- survfit(Surv(osm, ose) ~ 1,data=ilhd) 
fiten <- survfit(Surv(osm, ose) ~ 1,data=ilen) 
summary(ilen$osm)

nrm1 <-cuminc(ilhdos1y$aGVHD34m, ilhdos1y$agvhd34n, rho=0, cencode=0,
              na.action=na.omit)
timepoints(nrm1, times=100/30)



os.cox <- npsurv(cuminc(aGVHD34m , agvhd34n)~1,data=ilhdos1y)
ilhd$nrme
ilhd$nrmm
a<-etmCIF(Surv(aGVHD34m, agvhd34)~1, etype =agvhd34, data=ilhdos1y)
summary(a,times=c(100/30,12),level = 0.95)

a<-etmCIF(Surv(aGVHDm, GVHD)~1, etype =GVHD, data=ilhdos1y)
summary(a,times=c(100/30,12),level = 0.95)
ilhd$GVHD


a<-etmCIF(Surv(aGVHD34m, agvhd34n)~1, etype =agvhd34n, data=ilenos1y)
summary(a,times=c(100/30,12),level = 0.95)

a<-etmCIF(Surv(aGVHDm, GVHD)~1, etype =GVHD, data=ilenos1y)
summary(a,times=c(100/30,12),level = 0.95)
min(ilenos1y$aGVHD34m)
summary(nrm1)
plot(nrm1)
?cuminc
#######################################Cox############################
ilhd <- mutate(tbl_df(read.csv('IL18FINAL.CSV')),il18prelog=log2(freeIL18pre),age=Age,cr1=as.factor(cr),freeIL18d28=as.numeric(freeIL18d28),
               Gratwohl2vsrest=Gratwohl2vsrest==2, mismatch3rest1=mismatch3rest1==3, donorsex=donorsex==2,sex=sex==2,myeloid1lymphoid2=myeloid1lymphoid2==2,
               atg=atg==1,ric2rest1=ric2rest1==2,IL18BPaprelog=log2(IL18BPapre), allil18prelog=log2(IL18pre))

hdos1y.cox<- coxph(Surv(osm, ose)~il18prelog,data=ilhdos1y)
hdos1y.cox<- coxph(Surv(aGVHDm, GVHD)~il18prelog,data=ilhdos1y)
hdos1y.cox<- coxph(Surv(aGVHD34m, agvhd34)~il18prelog,data=ilhdos1y)
hdos1y.cox <- coxph(Surv(osm, ose)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                    +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilhdos1y)

ilhda<- filter(ilhd, GVHD==1)%>%mutate(IL18d28log=log2(IL18d28), freeIL18d28log=log2(freeIL18d28)
                                       ,IL18BPad28log=log2(IL18BPad28))
ilhdosasplit <- survSplit(Surv(osa, ose) ~ ., data= ilhda, cut=c(12,125),
                         episode= "tgroup", id="id")
ilhdosa1y <- filter(ilhdosasplit,tgroup==1)
hdos1y.cox<- coxph(Surv(osa, ose)~il18prelog,data=ilhdosa1y )
summary(hdos1y.cox)
ilhda<- filter(ilhd, GVHD==1,nrma>0)%>%mutate(IL18d28log=log2(IL18d28), freeIL18d28log=log2(freeIL18d28)
                                              ,IL18BPad28log=log2(IL18BPad28))
ilhdnrmasplit <- survSplit(Surv(nrma, nrme) ~ ., data= ilhda, cut=c(12,125),
                          episode= "tgroup", id="id")
ilhdnrma1y <- filter(ilhdnrmasplit,tgroup==1)
hdnrma1y.cox<- coxph(Surv(nrma, nrme)~ il18prelog,data=ilhdnrma1y)
ilhdnrmasplit <- survSplit(Surv(nrma, relapse) ~ ., data= ilhda, cut=c(12,125),
                           episode= "tgroup", id="id")
ilhdnrma1y <- filter(ilhdnrmasplit,tgroup==1)
hdnrma1y.cox<- coxph(Surv(nrma,relapse)~il18prelog,data=ilhdnrma1y)
summary(hdnrma1y.cox)
######################table#####################################
ilhd <- mutate(tbl_df(read.csv('IL18HD.CSV')),il18prelog=log2(freeIL18pre),age=Age,cr1=as.factor(cr),freeIL18d28=as.numeric(freeIL18d28),
               Gratwohl2vsrest=Gratwohl2vsrest==2, mismatch3rest1=mismatch3rest1==3, donorsex=donorsex==2,sex=sex==2,myeloid1lymphoid2=myeloid1lymphoid2==2,
               atg=atg==1,ric2rest1=ric2rest1==2)
ilhdossplit <- survSplit(Surv(osm, ose) ~ ., data= ilhd, cut=c(60,125),
                         episode= "tgroup", id="id")

ilhdos1y <- mutate(filter(ilhdossplit,tgroup==1),agvhd34=GVHDgrade>1,agvhd34n=agvhd34*1,aGVHD34m=ifelse(GVHDgrade>1,aGVHDm,osm))
hdnrmm1y.cox<- coxph(Surv(osm, ose)~il18prelog,data=ilhdos1y)
table(ilhdos1y$ric2rest1 ,ilhdos1y$ose)
summary(ilhdos1y$Gratwohl2vsrest)
ilhdnrmmsplit <- survSplit(Surv(nrmm, nrme) ~ ., data= ilhd, cut=c(60,100),
                         episode= "tgroup", id="id")
ilhdnrmm1y <-  mutate(filter(ilhdnrmmsplit,tgroup==1))
hdnrmm1y.cox<- coxph(Surv(nrmm, nrme)~il18prelog>8.5,data=ilhdnrmm1y)
hdnrmm1y.cox<- coxph(Surv(nrmm, relapse)~il18prelog,data=ilhdnrmm1y)

hdnrm.cox <- coxph(Surv(nrmm, nrme)~il18prelog+Age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                   +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilhdnrmm1y)
summary(hdnrm.cox)
hdrelapse.cox <- coxph(Surv(nrmm, relapse)~il18prelog+Age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                   +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilhdnrmm1y)
summary(hdrelapse.cox )

table(ilhdnrmm1y$ric2rest1,ilhdnrmm1y$nrme)
table(ilhdnrmm1y$ric2rest1,ilhdnrmm1y$relapse)

#################spearman#################
rcorr(ilhd$MIGpre, ilhd$IL18BPapre,type = 'spearman')

cor.test(ilhd$IL18BPapre, ilhd$IL18pre,
         method =  "spearman",
         conf.level = 0.95, continuity = FALSE)

rcorr(ilen$IL18BPapre, ilen$IL18BPapre,type = 'spearman')
cor.test(ilen$CXCL9pre, ilen$IL18BPapre,
         method =  "spearman",
         conf.level = 0.95, continuity = FALSE)

##############################ESSEN#################################################
ilen <- mutate(tbl_df(read.csv('IL18EN.CSV')),il18prelog=log2(freeIL18pre),cr1=as.factor(cr),
               Gratwohl2vsrest=Gratwohl2vsrest==2, mismatch3rest1=mismatch3rest1==3, donorsex=donorsex==2,
               sex=sex==2,myeloid1lymphoid2=myeloid1lymphoid2==2,IL18BPaprelog=log2(IL18BPapre), allil18prelog=log2(IL18pre),
               atg=atg==1,ric2rest1=ric2rest1==2,nrma=nrmm-aGVHDm,osan=osm-aGVHDm)
hist(ilen$il18prelog)
summary(ilen$CXCL9pre)

ilenossplit <- survSplit(Surv(osm, ose) ~ ., data= ilen, cut=c(60,125),
                         episode= "tgroup", id="id")
ilenos1y <- mutate(filter(ilenossplit,tgroup==1),agvhd34=GVHDgrade>1,aGVHD34m=ifelse(GVHDgrade>1,aGVHDm,osm),agvhd34n=agvhd34*1)
enos1y.cox<- coxph(Surv(osm, ose)~il18prelog,data=ilenos1y)
enos1y.cox<- coxph(Surv(aGVHDm, GVHD)~il18prelog,data=ilenos1y)
enos1y.cox<- coxph(Surv(aGVHD34m, agvhd34)~il18prelog,data=ilenos1y)
summary(enos1y.cox)

ilena<- filter(ilen, GVHD==1,osa>0)
ilenosasplit <- survSplit(Surv(osa, ose) ~ ., data= ilena, cut=c(12,125),
                          episode= "tgroup", id="id")
ilenosa1y <- filter(ilenosasplit,tgroup==1)
rnos1y.cox<- coxph(Surv(osa, ose)~il18prelog,data=ilenosa1y )
summary(rnos1y.cox)
ilenanrm<- filter(ilen, GVHD==1,nrma>0)
ilennrmasplit <- survSplit(Surv(nrma, nrme) ~ ., data= ilenanrm, cut=c(12,125),
                           episode= "tgroup", id="id")
ilennrma1y <- filter(ilennrmasplit,tgroup==1,nrma>0)
hdnrma1y.cox<- coxph(Surv(nrma, nrme)~il18prelog,data=ilennrma1y )
ilennrmasplit <- survSplit(Surv(nrma, relapse) ~ ., data= ilenanrm, cut=c(12,125),
                           episode= "tgroup", id="id")
ilennrma1y <- filter(ilennrmasplit,tgroup==1)
hdnrma1y.cox<- coxph(Surv(nrma,relapse)~il18prelog,data=ilennrma1y )
summary(hdnrma1y.cox)

enos.cox <- coxph(Surv(osm, ose)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                  +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilenos1y)
summary(enos.cox)

table(ilenos1y$ric2rest1,ilenos1y$ose)
table(ilennrmm1y$ric2rest1,ilennrmm1y$nrme)
table(ilennrmm1y$ric2rest1,ilennrmm1y$relapse)

summary(ilenos1y$ric2rest1)

ilennrmmsplit <- survSplit(Surv(nrmm, nrme) ~ ., data= ilen, cut=c(60,100),
                           episode= "tgroup", id="id")

ilennrmm1y <-  mutate(filter(ilennrmmsplit,tgroup==1))
hdnrmm1y.cox<- coxph(Surv(nrmm, nrme)~il18prelog,data=ilennrmm1y)
hdnrmm1y.cox<- coxph(Surv(nrmm, relapse)~il18prelog,data=ilennrmm1y)
summary(hdnrmm1y.cox)

enos.cox <- coxph(Surv(nrmm, nrme)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                  +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilennrmm1y)
enos.cox <- coxph(Surv(relapsem, relapse)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                  +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilennrmm1y)
summary(enos.cox)
#######################validation###############################
summary(ilhd$il18prelog)
ilhdossplit <- survSplit(Surv(osm, ose) ~ ., data= ilhd, cut=c(60,125),
                         episode= "tgroup", id="id")
ilhdos1y <- filter(ilhdossplit,tgroup==1)%>%mutate(il18precut=il18prelog>8.876,il18precut2=il18prelog> 9.457)
hdos1y.cox <- coxph(Surv(osm, ose)~il18prelog>9.67,data=ilhdos1y)
hdos1y.cox <- coxph(Surv(osm, ose)~il18precut2,data=ilhdos1y)
summary(hdos1y.cox )
hdos1y.SURV <- survdiff(Surv(osm, ose)~il18prelog>9.67,data=ilhdos1y)

hdos1y.SURV$chisq

mstat <- maxstat.test(Surv(osm, ose)~ il18prelog
                      , data=ilhdos1y,
                      smethod="LogRank",pmethod="exactGauss",
                      abseps=0.01)
plot(mstat)
abline(v=8.5)
hdos1y.cox <- coxph(Surv(osm, ose)~il18prelog2+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                    +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilhdos1y)

hdos1y.cox <- coxph(Surv(osm, ose)~il18precut2+Age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                    +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilhdos1y)

summary(hdos1y.cox)
cox.zph(hdos1y.cox)

ilenossplit <- survSplit(Surv(osm, ose) ~ ., data= ilen, cut=c(60,125),
                         episode= "tgroup", id="id")

ilenos1y <- mutate(filter(ilenossplit,tgroup==1), off1=il18prelog*0.220309,il18precut=il18prelog>8.876,il18precut2=il18prelog> 9.457)%>%mutate(off2=il18precut*0.296946,
                                                                                                                                               off3=il18precut2*0.363761)
enos1y.cox <- coxph(Surv(osm, ose)~il18precut2,data=ilenos1y )

enos.cox <- coxph(Surv(osm, ose)~il18prelog+offset(off1)+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                  +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilenos1y)
summary(enos1y.cox)
enos.cox <- coxph(Surv(osm, ose)~il18precut2+offset(off3)+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                  +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilenos1y)
summary(enos.cox)
###nrm********
ilhdnrmmsplit <- survSplit(Surv(nrmm, relapse) ~ ., data= ilhd, cut=c(60,125),
                           episode= "tgroup", id="id")
ilhdnrmm1y <-  filter(ilhdnrmmsplit,tgroup==1)%>%mutate(freeIL18d28log=log2(freeIL18d28),il18precut=il18prelog>8.876,il18precut2=il18prelog> 9.457)
hdnrmm1y.cox <- coxph(Surv(nrmm, nrme)~il18precut2,data=ilhdnrmm1y)
summary(hdnrmm1y.cox)
mstat <- maxstat.test(Surv(nrmm, nrme)~ il18prelog,
                      data=ilhdnrmm1y,
                      smethod="LogRank",
                      abseps=0.01)
plot(mstat)
abline(v=9.457)
summary(ilhd$freeIL18pre)
hdnrmm1y.cox <- coxph(Surv(nrmm, relapse)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                    +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilhdnrmm1y)
hdnrmm1y.cox <- coxph(Surv(nrmm, nrme)~il18precut+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                      +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilhdnrmm1y)
hdnrmm1y.cox <- coxph(Surv(nrmm, nrme)~il18precut2+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                      +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilhdnrmm1y)
summary(hdnrmm1y.cox)
cox.zph(hdnrmm1y.cox)
plot(cox.zph(hdnrmm1y.cox)[1])
abline(0,0, col=2,lwd=2)

ilennrmmsplit <- survSplit(Surv(nrmm, nrme) ~ ., data= ilen, cut=c(60,125),
                         episode= "tgroup", id="id")
ilennrmm1y <- mutate(filter(ilennrmmsplit,tgroup==1), off1=il18prelog*0.33055,il18precut=il18prelog>8.876,il18precut2=il18prelog> 9.457)%>%mutate(off2=il18precut*0.74191,          
                                                                                                                                                  off3=il18precut2*0.57543)

ennrmm1y.cox <- coxph(Surv(nrmm, nrme)~il18precut2,data=ilennrmm1y)
summary(ennrmm1y.cox)
mstaten <- maxstat.test(Surv(nrmm, nrme)~ il18prelog
                      , data=ilennrmm1y,
                      smethod="LogRank",
                      abseps=0.01)
plot(mstaten)
abline(v=8.695)
summary(ilen$il18prelog)
summary(ilen$freeIL18pre)
ennrmm1y.cox <- coxph(Surv(nrmm, nrme)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                  +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilennrmm1y)

ennrmm1y.cox <- coxph(Surv(nrmm, nrme)~il18prelog+offset(off1) +age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                      +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilennrmm1y)
ennrmm1y.cox <- coxph(Surv(nrmm, nrme)~il18precut+offset(off2) +age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                      +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilennrmm1y)
ennrmm1y.cox <- coxph(Surv(nrmm, nrme)~il18precut2+offset(off3) +age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                      +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilennrmm1y)
summary(ennrmm1y.cox)
ennrmm1y.cox
cox.zph(ennrmm1y.cox)
ennrmm1y.cox <- coxph(Surv(nrmm, nrme)~il18prelog+offset(off1)+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                  +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilennrmm1y)
summary(ennrmm1y.cox)
############PE validation#################################

ilennrmm1y <- mutate(filter(ilennrmmsplit,tgroup==1), off1=il18prelog*(0.33055))
ModelsNRM <- list ("Essen"=coxph(Surv(nrmm, nrme)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                                +sex+myeloid1lymphoid2+atg+ric2rest1, data=ilennrmm1y,x=TRUE),
                  "offset"=coxph(Surv(nrmm, nrme)~offset(off1)+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                                +sex+myeloid1lymphoid2+atg+ric2rest1, data=ilennrmm1y,x=TRUE))


PredError.632plusleee<- pec(object=ModelsNRM ,
                            formula=Surv(nrmm, nrme)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                            +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilennrmm1y, reference = FALSE,
                            exact=TRUE,
                            cens.model="cox",
                            splitMethod="Boot632plus",
                            B=100,
                            verbose=TRUE)
plot(PredError.632plusleee)

plot(PredError.632plusleee,lty = c(1,2),col = c('red','black'), xlim=c(0,59),xlab = 'Time After alloSCT(months)',legend = TRUE)
legend('bottomright', c("loguni.ref","loguni",'logmulti.ref', "logmuti") , 
       lty = c(1,2,1,2), col=c('black','green','blue','red'), bty='n', cex=0.9,lwd = 2, y.intersp=2)

ilenos1y <- mutate(filter(ilenossplit,tgroup==1), off1=il18prelog*0.220309)

Modelsos <- list ("Essen"=coxph(Surv(osm,ose)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                                 +sex+myeloid1lymphoid2+atg+ric2rest1, data=ilenos1y,x=TRUE),
                   "offset"=coxph(Surv(osm,ose)~offset(off1)+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                                  +sex+myeloid1lymphoid2+atg+ric2rest1, data=ilenos1y,x=TRUE))

PredError.632plusos<- pec(object=Modelsos ,
                            formula=Surv(osm,ose)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                            +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilenos1y, reference = FALSE,
                            exact=TRUE,
                            cens.model="cox",
                            splitMethod="Boot632plus",
                            B=100,
                            verbose=TRUE)
plot(PredError.632plusos,lty = c(1,2),col = c('red','black'), xlim=c(0,59),xlab = 'Time After alloSCT(months)',legend = TRUE) 


ilhd <- mutate(tbl_df(read.csv('IL18HD.CSV')),il18prelog=log2(freeIL18pre),il18d28log=log2(freeIL18d28),age=Age,cr1=as.factor(cr),freeIL18d28=as.numeric(freeIL18d28))
ilhdnrmmsplit <- survSplit(Surv(nrmm, nrme) ~ ., data= ilhd, cut=c(60,100),
                           episode= "tgroup", id="id")
ilhdnrmm1y <-  filter(ilhdnrmmsplit,tgroup==1)
ilhdnrmm1y
ModelsNRMhd <- list ("HD"=coxph(Surv(nrmm, nrme)~il18prelog, data=ilhdnrmm1y,x=TRUE),
                     "multi"=coxph(Surv(nrmm, nrme)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                                           +sex+myeloid1lymphoid2+atg+ric2rest1, data=ilhdnrmm1y,x=TRUE),
                     "multiref"=coxph(Surv(nrmm, nrme)~age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                                    +sex+myeloid1lymphoid2+atg+ric2rest1, data=ilhdnrmm1y,x=TRUE))
ApparrentCindexmodelch <- pec::cindex(ModelsNRMhd,
                                      formula=Surv(nrmm, nrme)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                                      +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilhdnrmm1y,
                                      eval.times=seq(0,60,5))
plot(ApparrentCindexmodelch,lty = c(1,1,1),col=c('red','brown','blue') , xlim=c(0,60),xlab = 'Time After alloSCT(months)',legend =TRUE)


ilhdossplit <- survSplit(Surv(osm, ose) ~ ., data= ilhd, cut=c(60,125),
                         episode= "tgroup", id="id")

ilhdos1y <- mutate(filter(ilhdossplit,tgroup==1),agvhd34=GVHDgrade>1,agvhd34n=agvhd34*1,aGVHD34m=ifelse(GVHDgrade>1,aGVHDm,osm))

ModelsOSMhd <- list ("HD"=coxph(Surv(osm, ose)~il18prelog, data=ilhdos1y,x=TRUE),
                     "multi"=coxph(Surv(osm, ose)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                                   +sex+myeloid1lymphoid2+atg+ric2rest1, data=ilhdos1y,x=TRUE),
                     "multiref"=coxph(Surv(osm, ose)~age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                                      +sex+myeloid1lymphoid2+atg+ric2rest1, data=ilhdos1y,x=TRUE))
PredErrorhdos<- pec(object=ModelsOSMhd ,
                            formula=Surv(osm, ose)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                    +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilhdos1y, reference = TRUE,
                            exact=TRUE,
                            cens.model="cox",
                            splitMethod="non",
                            B=100,
                            verbose=TRUE)
ApparrentCindexmodelch <- pec::cindex(ModelsOSMhd,
                                      formula=Surv(osm, ose)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                                      +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilhdos1y,
                                      eval.times=seq(0,60,5))
plot(ApparrentCindexmodelch,lty = c(1,1,1),col=c('red','brown','blue') , xlim=c(0,60),xlab = 'Time After alloSCT(months)',legend =TRUE)



ilen <- mutate(tbl_df(read.csv('IL18EN.CSV')),il18prelog=log2(freeIL18pre),cr1=as.factor(cr),
               Gratwohl2vsrest=Gratwohl2vsrest==2, mismatch3rest1=mismatch3rest1==3, donorsex=donorsex==2,
               sex=sex==2,myeloid1lymphoid2=myeloid1lymphoid2==2,IL18BPaprelog=log2(IL18BPapre), allil18prelog=log2(IL18pre),
               atg=atg==1,ric2rest1=ric2rest1==2)

ilenossplit <- survSplit(Surv(osm, ose) ~ ., data= ilen, cut=c(60,125),
                         episode= "tgroup", id="id")
ilenos1y <- mutate(filter(ilenossplit,tgroup==1),agvhd34=GVHDgrade>1,aGVHD34m=ifelse(GVHDgrade>1,aGVHDm,osm),agvhd34n=agvhd34*1)%>%filter(aGVHD34m>0,aGVHDm>0)

Modelsosmen <- list ("Essen"=coxph(Surv(osm, ose)~il18prelog, data=ilenos1y,x=TRUE),
                     "multi"=coxph(Surv(osm, ose)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                                   +sex+myeloid1lymphoid2+atg+ric2rest1, data=ilenos1y,x=TRUE),
                     "multiref"=coxph(Surv(osm, ose)~age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                                      +sex+myeloid1lymphoid2+atg+ric2rest1, data=ilenos1y,x=TRUE))
PredErrorenos<- pec(object=Modelsosmen ,
                    formula=Surv(osm, ose)~1,data=ilenos1y, reference = TRUE,
                    exact=TRUE,
                    cens.model="cox",
                    splitMethod="non",
                    B=100,
                    verbose=TRUE)
plot(PredErrorenos,xlim = c(0,59))

ApparrentCindexmodelch <- pec::cindex(Modelsosmen,
                                      formula=Surv(osm, ose)~1,data=ilenos1y,
                                      eval.times=seq(0,60,5))
plot(ApparrentCindexmodelch ,lty = c(1,1,1),col=c('red','brown','blue') , xlim=c(0,60),xlab = 'Time After alloSCT(months)',legend =TRUE)


ilennrmmsplit <- survSplit(Surv(nrmm, nrme) ~ ., data= ilen, cut=c(60,125),
                         episode= "tgroup", id="id")
ilennrmm1y <- mutate(filter(ilennrmmsplit,tgroup==1),agvhd34=GVHDgrade>1,aGVHD34m=ifelse(GVHDgrade>1,aGVHDm,osm),agvhd34n=agvhd34*1)%>%filter(aGVHD34m>0,aGVHDm>0)

ModelsNRMen <- list ("HD"=coxph(Surv(nrmm, nrme)~il18prelog, data=ilennrmm1y,x=TRUE),
                     "multi"=coxph(Surv(nrmm, nrme)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                                   +sex+myeloid1lymphoid2+atg+ric2rest1, data=ilennrmm1y,x=TRUE),
                     "multiref"=coxph(Surv(nrmm, nrme)~age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                                      +sex+myeloid1lymphoid2+atg+ric2rest1, data=ilennrmm1y,x=TRUE))
PredErrorennrm<- pec(object=ModelsNRMen,
                    formula=Surv(nrmm, nrme)~1,data=ilennrmm1y, reference = TRUE,
                    exact=TRUE,
                    cens.model="cox",
                    splitMethod="non",
                    B=100,
                    verbose=TRUE)
plot(PredErrorennrm,xlim = c(0,59))

ApparrentCindexmodelch2<- pec::cindex(ModelsNRMen,
                                      formula=Surv(nrmm, nrme)~1,data=ilennrmm1y,
                                      eval.times=seq(0,60,5))
plot(ApparrentCindexmodelch2,lty = c(1,1,1),col=c('red','brown','blue') , xlim=c(0,60),xlab = 'Time After alloSCT(months)',legend =TRUE)
 ##############IL18PRE VS IL18 DAY 28################################   

ilhdnrmmsplit <- survSplit(Surv(osm,ose) ~ ., data= ilhd, cut=c(12,100),
                           episode= "tgroup", id="id")
ilhdnrmm1y <-  filter(ilhdnrmmsplit,tgroup==1)

ilhdnrmm1ypepred28 <-  filter(ilhdnrmm1y,!is.na(freeIL18pre),!is.na(freeIL18d28))%>%mutate(il18prelog=log2(freeIL18pre),
                                                                                           il18d28log=log2(freeIL18d28))

ModelsNRMhd <- list ("IL18pre"=coxph(Surv(nrmm, nrme)~il18prelog, data=ilhdnrmm1ypepred28,x=TRUE),
                     "IL18D28"=coxph(Surv(nrmm, nrme)~il18d28log, data=ilhdnrmm1ypepred28,x=TRUE))

ModelsNRMhd <- list ("IL18pre"=coxph(Surv(osm, ose)~il18prelog, data=ilhdnrmm1ypepred28,x=TRUE),
                     "IL18D28"=coxph(Surv(osm, ose)~il18d28log, data=ilhdnrmm1ypepred28,x=TRUE))

PredError.632plusleee1<- pec(object=ModelsNRMhd ,
                            formula=Surv(nrmm, nrme)~1,data=ilhdnrmm1ypepred28, reference = FALSE,
                            exact=TRUE,
                            cens.model="marginal",
                            splitMethod="non",
                            B=100,
                            verbose=TRUE)
plot(PredError.632plusleee1,xlim=c(0,11.9))
ilhdnrmm1y.na<- filter(ilhdnrmm1y, !is.na(nrmm),!is.na(age),!is.na(Gratwohl2vsrest),!is.na(mismatch3rest1),!is.na(donorsex),
                       !is.na(sex),!is.na(myeloid1lymphoid2),!is.na(atg),!is.na(ric2rest1))
ApparrentCindexmodelch <- pec::cindex(ModelsNRMhd,
                                      formula=Surv(nrmm, nrme)~il18prelog+age+Gratwohl2vsrest+mismatch3rest1+donorsex 
                                      +sex+myeloid1lymphoid2+atg+ric2rest1,data=ilhdnrmm1y,
                                      eval.times=seq(0,60,5))
plot(ApparrentCindexmodelch,lty = c(1,2),col=c('red','red','black') , xlim=c(0,60),xlab = 'Time After alloSCT(months)',legend =TRUE)


##################PLOTS############################

ilhd$Qfreeil18pre<-CutQ(ilhd$il18prelog, breaks=quantile(ilhd$il18prelog, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)


ilhd$con14
ilhdcut<- mutate(ilhd, il18precut=il18prelog>8.876,il18precut2=il18prelog> 9.457)
fitsurvnrm <- npsurv(Surv(nrmm, factor(con14))~Qfreeil18pre,data=ilhd)

fitsurvnrm <- npsurv(Surv(nrmm, factor(cr))~il18precut2,data=ilhdcut)
par(mgp=c(1.5, 1, 0))
survplot(fitsurvnrm ,state=2,n.risk=TRUE, xlim=c(0, 60),ylim=c(0, 1),conf.int=0.1,time.inc=10,label.curves=TRUE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('red','black'),
         lty = c(1,2,2,1),xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)
legend('topright', c("Q1","Q2","Q3","Q4") , 
       lty = c(1,2,2,1), col=c('red','black'), bty='n', cex=0.9,lwd = 2, y.intersp=2)

fitsurvos <- npsurv(Surv(osm,ose)~Q2,data=ilhd )
fitsurvos <- npsurv(Surv(osm,ose)~il18precut2,data=ilhdcut)
par(mgp=c(1.5, 1, 0))
survplot(fitsurvos  ,n.risk=TRUE, xlim=c(0, 60),ylim=c(0, 1),conf.int=0.1,time.inc=10,label.curves=TRUE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('red','black'),
         lty = c(1,2,2,1),xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)
legend('topright', c("Q1","Q2","Q3","Q4") , 
       lty = c(1,2,2,1), col=c('red','black'), bty='n', cex=0.9,lwd = 2, y.intersp=2)

ilhda<- filter(ilhd, GVHD==1)
ilhda$Qfreeil18pre <- CutQ(ilhda $il18prelog, breaks = quantile(ilhda$il18prelog, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhda$Qil18pre <- CutQ(ilhda $IL18pre, breaks = quantile(ilhda$IL18pre, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhda$Qil18pa <- CutQ(ilhda$IL18BPapre, breaks = quantile(ilhda$IL18BPapre, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhda$QIL18d28 <- CutQ(ilhda$IL18d28, breaks = quantile(ilhda$IL18d28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhda$QfreeIL18d28 <- CutQ(ilhda$freeIL18d28, breaks = quantile(ilhda$freeIL18d28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhda$QIL18BPad28 <- CutQ(ilhda$IL18BPad28, breaks = quantile(ilhda$IL18BPad28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)

fitsurvnrm <- npsurv(Surv(nrma, factor(cra))~QfreeIL18d28,data=ilhda)
par(mgp=c(1.5, 1, 0))
survplot(fitsurvnrm ,state=1,n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=3,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('red','black'),
         lty = c(1,2,2,1),xlab='Time After aGVHD (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)
legend('topright', c("Q1","Q2","Q3","Q4") , 
       lty = c(1,2,2,1), col=c('red','black'), bty='n', cex=0.9,lwd = 2, y.intersp=2)

fitsurvosma <- npsurv(Surv(osa,ose)~QfreeIL18d28,data=ilhda)
par(mgp=c(1.5, 1, 0))
survplot(fitsurvosma,n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=3,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('red','black'),
         lty = c(1,2,2,1),xlab='Time After aGVHD (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)

ilhd <- mutate(tbl_df(read.csv('IL18FINAL.CSV')),il18prelog=log2(freeIL18pre),age=Age,cr1=as.factor(cr),freeIL18d28=as.numeric(freeIL18d28),
               Gratwohl2vsrest=Gratwohl2vsrest==2, mismatch3rest1=mismatch3rest1==3, donorsex=donorsex==2,sex=sex==2,myeloid1lymphoid2=myeloid1lymphoid2==2,
               atg=atg==1,ric2rest1=ric2rest1==2,IL18BPaprelog=log2(IL18BPapre), allil18prelog=log2(IL18pre))
ilhdnrmmsplit <- survSplit(Surv(nrmm, nrme) ~ ., data= ilhd, cut=c(60,125),
                           episode= "tgroup", id="id")
ilhdnrmm1y <- mutate(filter(ilhdnrmmsplit,tgroup==1),agvhd34=GVHDgrade>1,aGVHD34m=ifelse(GVHDgrade>1,aGVHDm,osm),agvhd34n=agvhd34*1)%>%filter(aGVHD34m>0,aGVHDm>0)

ilhdnrmmd28<- mutate(ilhd,mrmm28=nrmm-(28/30),osm28=osm-(28/30),IL18d28log=log2(IL18d28), freeIL18d28log=log2(freeIL18d28),
IL18BPad28log=log2(IL18BPad28)) %>%filter(mrmm28>=0,!is.na(freeIL18pre))

ilhdnrmmsplit28 <- survSplit(Surv(mrmm28,relapse) ~ ., data= ilhdnrmmd28, cut=c(12,125),
                           episode= "tgroup", id="id1")
ilhdnrmm1y28 <- mutate(filter(ilhdnrmmsplit28,tgroup==1,!is.na(freeIL18pre)),agvhd34=GVHDgrade>1,aGVHD34m=ifelse(GVHDgrade>1,aGVHDm,osm),agvhd34n=agvhd34*1)%>%filter(aGVHD34m>0,aGVHDm>0)


hdnrmmd28.cox<- coxph(Surv(mrmm28, nrme)~freeIL18d28log,data=ilhdnrmm1y28)
hdnrmmd28.cox<- coxph(Surv(mrmm28, relapse)~freeIL18d28log,data=ilhdnrmmd28)
summary(hdnrmmd28.cox)

ilhdosmsplit <- survSplit(Surv(osm, ose) ~ ., data= ilhd, cut=c((12+(28/30)),125),
                           episode= "tgroup", id="id")
ilhdos1y <- mutate(filter(ilhdosmsplit,tgroup==1),agvhd34=GVHDgrade>1,aGVHD34m=ifelse(GVHDgrade>1,aGVHDm,osm),agvhd34n=agvhd34*1)%>%filter(aGVHD34m>0,aGVHDm>0)

ilhdnrmmd28$QIL18d28 <- CutQ(ilhdnrmmd28$IL18d28, breaks = quantile(ilhdnrmmd28$IL18d28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhdnrmmd28$QfreeIL18d28 <- CutQ(ilhdnrmmd28$freeIL18d28, breaks = quantile(ilhdnrmmd28$freeIL18d28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhdnrmmd28$QIL18BPad28 <- CutQ(ilhdnrmmd28$IL18BPad28, breaks = quantile(ilhdnrmmd28$IL18BPad28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhdnrmmd28$mrmm28
fitsurvosma <- npsurv(Surv(mrmm28,factor(cr))~QfreeIL18d28,data=ilhdnrmmd28)
par(mgp=c(1.5, 1, 0))
survplot(fitsurvosma,state=1,n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('red','black'),
         lty = c(1,2,2,1),xlab='Time After day28 (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)

ilhdosmd28<- mutate(ilhd,mrmm28=nrmm-(28/30),osm28=osm-(28/30),IL18d28log=log2(IL18d28), freeIL18d28log=log2(freeIL18d28),
                    IL18BPad28log=log2(IL18BPad28))%>%filter(osm28>0,!is.na(freeIL18pre))
ilhdosmsplit <- survSplit(Surv(osm28, ose) ~ ., data= ilhdosmd28, cut=c(12,125),
                          episode= "tgroup", id="id1")
ilhdos1y <- mutate(filter(ilhdosmsplit,tgroup==1),agvhd34=GVHDgrade>1,aGVHD34m=ifelse(GVHDgrade>1,aGVHDm,osm),agvhd34n=agvhd34*1)%>%filter(aGVHD34m>0,aGVHDm>0)

hdosmd28.cox<- coxph(Surv(osm28, ose)~freeIL18d28log,data=ilhdos1y )
summary(hdosmd28.cox)

ilhdosmd28$QIL18d28 <- CutQ(ilhdosmd28$IL18d28, breaks = quantile(ilhdosmd28$IL18d28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhdosmd28$QfreeIL18d28 <- CutQ(ilhdosmd28$freeIL18d28, breaks = quantile(ilhdosmd28$freeIL18d28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhdosmd28$QIL18BPad28 <- CutQ(ilhdosmd28$IL18BPad28, breaks = quantile(ilhdosmd28$IL18BPad28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
fitsurvosma <- npsurv(Surv(osm28,ose)~QfreeIL18d28,data=ilhdosmd28)
par(mgp=c(1.5, 1, 0))
survplot(fitsurvosma,n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=2,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('red','black'),
         lty = c(1,2,2,1),xlab='Time After day28 (Months)', ylab='Survival Probability',cex.xlab=0.9)

ilhdagvhd34m<- mutate(ilhdosmd28,agvhd34=ifelse(GVHDgrade>1,'2',GVHDgrade),agvhdre=((GVHDgrade==3)*1),agvhd34csc=((GVHDgrade>1)*1),
                    aGVHD34m=ifelse(GVHDgrade>1,aGVHDm-(28/30),osm28), aGVHDrem=ifelse(GVHDgrade==3,aGVHDm-(28/30),osm28),IL18d28log=log2(IL18d28),
                    freeIL18d28log=log2(freeIL18d28),IL18BPad28log=log2(IL18BPad28))%>%filter(aGVHD34m>0)

ilhdagvhd34m$QIL18d28 <- CutQ(ilhdagvhd34m$IL18d28, breaks = quantile(ilhdagvhd34m$IL18d28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhdagvhd34m$QfreeIL18d28 <- CutQ(ilhdagvhd34m$freeIL18d28, breaks = quantile(ilhdagvhd34m$freeIL18d28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhdagvhd34m$QIL18BPad28 <- CutQ(ilhdagvhd34m$IL18BPad28, breaks = quantile(ilhdagvhd34m$IL18BPad28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)

ilhdagvhdm$agvhdre
fitsurvagvhdm <- npsurv(Surv(aGVHD34m,factor(agvhd34csc))~QIL18d28,data=ilhdagvhd34m)
par(mgp=c(1.5, 1, 0))
survplot(fitsurvagvhdm,state=1,n.risk=TRUE, xlim=c(0, 25),ylim=c(0, 1),conf.int=0.1,time.inc=5,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('red','black'),
         lty = c(1,2,2,1),xlab='Time After day28 (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)
fitcoxagvhd34m <- coxph(Surv(aGVHD34m,agvhd34csc)~IL18d28log,data=ilhdagvhd34m)
summary(fitcoxagvhd34m)

ilhdagvhdrem<- mutate(ilhdosmd28,agvhd34=ifelse(GVHDgrade>1,'2',GVHDgrade),agvhdre=((GVHDgrade==3)*1),agvhd34csc=((GVHDgrade>1)*1),
                      aGVHD34m=ifelse(GVHDgrade>1,aGVHDm-(28/30),osm28), aGVHDrem=ifelse(GVHDgrade==3,aGVHDm-(28/30),osm28),IL18d28log=log2(IL18d28),
                      freeIL18d28log=log2(freeIL18d28),IL18BPad28log=log2(IL18BPad28))%>%filter(aGVHDrem>0)

ilhdagvhdrem$QIL18d28 <- CutQ(ilhdagvhdrem$IL18d28, breaks = quantile(ilhdagvhdrem$IL18d28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhdagvhdrem$QfreeIL18d28 <- CutQ(ilhdagvhdrem$freeIL18d28, breaks = quantile(ilhdagvhdrem$freeIL18d28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhdagvhdrem$QIL18BPad28 <- CutQ(ilhdagvhdrem$IL18BPad28, breaks = quantile(ilhdagvhdrem$IL18BPad28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
fitsurvagvhdrem <- npsurv(Surv(aGVHDrem,factor(agvhdre))~QfreeIL18d28,data=ilhdagvhdrem)
par(mgp=c(1.5, 1, 0))
survplot(fitsurvagvhdrem,state=1,n.risk=TRUE, xlim=c(0, 25),ylim=c(0, 1),conf.int=0.1,time.inc=5,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('red','black'),
         lty = c(1,2,2,1),xlab='Time After day28 (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)
hdnrmmd28.cox<- coxph(Surv(aGVHDrem, agvhdre)~freeIL18d28log,data=ilhdagvhdrem)
summary(hdnrmmd28.cox)
ilhdtma<-filter(ilhd, tmam>(28/30))%>%mutate(tmam28=tmam-(28/30),IL18d28log=log2(IL18d28),
                                             freeIL18d28log=log2(freeIL18d28),IL18BPad28log=log2(IL18BPad28))
summary(ilhdtma$QfreeIL18d28)
ilhdtma$QfreeIL18d28 <- CutQ(ilhdtma$freeIL18d28, breaks = quantile(ilhdtma$freeIL18d28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilhdtma$QIL18d28 <- CutQ(ilhdtma$IL18d28, breaks = quantile(ilhdtma$IL18d28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
fitsurvtma <- npsurv(Surv(tmam,factor(tma))~QIL18d28,data=ilhdtma)
par(mgp=c(1.5, 1, 0))
survplot(fitsurvtma,state=1,n.risk=TRUE, xlim=c(0, 25),ylim=c(0, 1),conf.int=0.1,time.inc=5,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('red','black'),
         lty = c(1,2,2,1),xlab='Time After day28 (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)

hdnrmmd28.cox<- coxph(Surv(tmam, tma)~IL18d28log,data=ilhdtma)
summary(hdnrmmd28.cox)

ilhd$Qil18post <- CutQ(ilhd$freeIL18d28, breaks = quantile(ilhd$freeIL18d28, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)

fitsurvnrm <- npsurv(Surv(nrmm, factor(cr))~Qil18post,data=ilhd )
par(mgp=c(1.5, 1, 0))
survplot(fitsurvnrm ,state=2,n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=3,label.curves=TRUE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('red','black'),
         lty = c(1,2,2,1),xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)
legend('topright', c("Q1","Q2","Q3","Q4") , 
       lty = c(1,2,2,1), col=c('red','black'), bty='n', cex=0.9,lwd = 2, y.intersp=2)

hdnrmmd28.cox<- coxph(Surv(aGVHD34m, agvhd34csc)~freeIL18d28log,data=ilhdtma)
summary(hdnrmmd28.cox)

ilen$Q2<- CutQ(ilen$il18prelog, breaks = quantile(ilen$il18prelog, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilencut<- mutate(ilen,il18precut=il18prelog>8.876,il18precut2=il18prelog> 9.457)
fitsurvnrmen <- npsurv(Surv(nrmm, cr1)~il18precut2,data=ilencut)
fitsurvnrmen <- npsurv(Surv(nrmm, cr1)~Q2,data=ilen)
par(mgp=c(1.5, 1, 0))
survplot(fitsurvnrmen  ,state=2,n.risk=TRUE, xlim=c(0, 60),ylim=c(0, 1),conf.int=0.1,time.inc=10,label.curves=TRUE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('red','black'),
         lty = c(1,2,2,1),xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)
legend('topright', c("Q1","Q2","Q3","Q4") , 
       lty = c(1,2,2,1), col=c('red','black'), bty='n', cex=0.9,lwd = 2, y.intersp=2)

fitsurvosen <- npsurv(Surv(osm,ose)~Q2,data=ilen )
fitsurvosen <- npsurv(Surv(osm,ose)~il18precut2,data=ilencut)
par(mgp=c(1.5, 1, 0))
survplot(fitsurvosen  ,n.risk=TRUE, xlim=c(0, 60),ylim=c(0, 1),conf.int=0.1,time.inc=10,label.curves=TRUE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('red','black'),
         lty = c(1,2,2,1),xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)
legend('topright', c("Q1","Q2","Q3","Q4") , 
       lty = c(1,2,2,1), col=c('red','black'), bty='n', cex=0.9,lwd = 2, y.intersp=2)

ilena<- filter(ilen, GVHD==1,nrma>=0)

ilena$Qfreeil18pre <- CutQ(ilena $freeIL18pre, breaks = quantile(ilena$freeIL18pre, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
ilena$Qil18pre <- CutQ(ilena $IL18pre, breaks = quantile(ilena$IL18pre, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)

ilena$Qil18pa <- CutQ(ilena$IL18BPapre, breaks = quantile(ilena$IL18BPapre, seq(0, 1, by = 0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE)
summary(ilena$IL18BPapre)
fitsurvnrmaen <- npsurv(Surv(nrma, factor(cra))~Qil18pa,data=ilena)
par(mgp=c(1.5, 1, 0))
survplot(fitsurvnrmaen ,state=2,n.risk=TRUE, xlim=c(0, 12),ylim=c(0, 1),conf.int=0.1,time.inc=3,label.curves=FALSE,lwd = 2,
         y.n.risk=-0.5,cex.n.risk=0.7,col.fill='white',sep.n.risk=0.06, cex.axis=2,
         col=c('red','black'),
         lty = c(1,2,2,1),xlab='Time After AlloSCT (Months)', ylab='Cumulative Incidence',cex.xlab=0.9)
legend('topright', c("Q1","Q2","Q3","Q4") , 
       lty = c(1,2,2,1), col=c('red','black'), bty='n', cex=0.9,lwd = 2, y.intersp=2)

###########################flow plots#################
ilhdx<- filter(ilhd, !is.na(ilhd $freeIL18pre),!is.na(ilhd$freeIL18d28))
library(ggalluvial)
ilhdx$Qpre <- ordered(CutQ(ilhdx $freeIL18pre, breaks = quantile(ilhdx$freeIL18pre, seq(1, 0, by = -0.25), na.rm = TRUE),  labels = NULL, na.rm = FALSE),
levels = c("Q4", "Q3", "Q2",'Q1'))
ilhdx$Qpost <- CutQ(ilhdx $freeIL18d28 , breaks = quantile(ilhdx$freeIL18d28, seq(1, 0, by = -0.25), na.rm = TRUE),  labels =NULL, na.rm = FALSE)
ilhdAlluvia <- mutate(select(ilhdx,Qpre, Qpost),count=1)

summary(ilhdx$Qpost)
shapiro.test(ilhdx$freeIL18pre)

is_alluvia_form(as.data.frame(ilhdAlluvia), axes = 1:3, silent = TRUE)
?quantile
ggplot(as.data.frame(ilhdAlluvia),
       aes(axis1 = Qpre, axis2 = Qpost)) +
  geom_alluvium(aes(fill = Qpre), width = 1/12) +
  geom_stratum(width = 1/12, fill = 'black', color = "grey") +
  geom_label(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c("freeIL18pre", "freeIL18d28"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1")


library(ggsci)
scale_fill_palname()
aggv<-as.data.frame(aggregate(ilhdAlluvia$count,by=list(ilhdAlluvia$Qpre,ilhdAlluvia$Qpost),sum))
ggplot(as.data.frame(aggv),
       aes(y = x, axis1 = Group.1, axis2 = Group.2)) +
  geom_alluvium(aes(fill = Group.1), width = 1/12) +
  geom_stratum(width = 1/12, fill = 'black', color = "grey") +
  geom_label(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c("freeIL18pre", "freeIL18d28"), expand = c(.05, .05)) +scale_fill_jco()

titanic_wide <- data.frame(Titanic)
ggplot(data = titanic_wide,
       aes(axis1 = Class, axis2 = Sex, axis3 = Age,
           y = Freq)) +
  scale_x_discrete(limits = c("Class", "Sex", "Age"), expand = c(.1, .05)) +
  xlab("Demographic") +
  geom_alluvium(aes(fill = Survived)) +
  geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE) +
  theme_minimal() +
  ggtitle("passengers on the maiden voyage of the Titanic",
          "stratified by demographics and survival")

  
summary(ilhdx$freeIL18pre)
sum(aggv$x)
xtabs(formula=count~Qpre + Qpost, data=ilhdAlluvia)
aggv<-as.data.frame(aggregate(ilhdAlluvia$count,by=list(ilhdAlluvia$Qpre,ilhdAlluvia$Qpost),sum))
aggv<- mutate(aggv, weight=4*x/395)
library(googleVis)
s=gvisSankey(aggv, from = "Group.1", to = "Group.2",weight="weight")
x=plot(s)

library(export)
filen <- paste("YOUR_DIR/ggplot")
graph2eps(x=x,file=filen, aspectr=2, font = "Times New Roman", 
          height = 5, bg = "white")
plot(1,2,3)
setEPS()
postscript("whatever.eps")
plot(rnorm(100), main="Hey Some Data")
dev.off()
###############ROC############################
ilhdnrme1<- mutate(ilhd,nrm1e=ifelse(nrmm<=12,nrme,0), os1e=ifelse(osm<=12,ose,0),il18prelog=log2(freeIL18pre),freeIL18d28log=log2(freeIL18d28))
ilhdnrme1d28<-mutate(ilhd,nrm1e=ifelse(nrmm<=12,nrme,0), nrm1ed28=ifelse(nrmm<=(12+(28/30.5)),nrme,0),os1ed28=ifelse(osm<=(12+(28/30.5)),ose,0),
                     il18prelog=log2(freeIL18pre),freeIL18d28log=log2(freeIL18d28))%>%filter(osm>=28/30.5)
ilhdnrmrocpre <-roc(ilhdnrme1$os1e, ilhdnrme1$il18prelog,na.rm=TRUE )
ilhdnrmrocd28 <-roc(ilhdnrme1d28$os1ed28, ilhdnrme1d28$freeIL18d28log,na.rm=TRUE )
plot(ilhdnrmrocpre)
lines(ilhdnrmrocd28, type="l",  col="red", bg="grey")
roc.test(ilhdnrmrocpre,ilhdnrmrocd28)
##################maxstat######################
library(cr17)
ilhd <- mutate(tbl_df(read.csv('IL18FINAL.CSV')),il18prelog=log2(freeIL18pre),age=Age,cr1=as.factor(cr),freeIL18d28=as.numeric(freeIL18d28),
               Gratwohl2vsrest=Gratwohl2vsrest==2, mismatch3rest1=mismatch3rest1==3, donorsex=donorsex==2,sex=sex==2,myeloid1lymphoid2=myeloid1lymphoid2==2,
               atg=atg==1,ric2rest1=ric2rest1==2,IL18BPaprelog=log2(IL18BPapre), allil18prelog=log2(IL18pre),il18cut=il18prelog>8.876118)
maxnrm<-mutate(ilhd, nrmmn=ifelse(nrmm>60,60,nrmm), nrmen=ifelse(nrmm>60&nrme==1,0,nrme),crn=ifelse(nrmm>60&cr>1,0,cr))%>%filter(!is.na(nrmmn),!is.na(crn),!is.na(il18prelog))
coxmax<- coxph(Surv(nrmmn,nrmen)~il18prelog,data=maxnrm)

x<-maxnrm[order(maxnrm$il18prelog),]
a<-x$il18prelog
length(maxnrm$crn)
p<-c(1:length(a))
q<-c(1:length(a))

median(maxnrm$il18prelog)

for( i in (2:588)){
  q[i]<-fitCuminc(time = x$nrmmn, risk = x$crn, group = (a-a[i])>0, cens = 0)$Tests$stat[2]
}
plot(a,q,ylim = c(0,12),xlim = c(8,10),ylab = c('Gray's test statistic'),xlab = c('log2(pre-transplant free IL-18)'))

plot(a,q,ylim = c(0,13),xlim = c(8,10),ylab = c('ln(Gray's test statistic)'),xlab = c('log2(pre-transplant free IL-18)'))

lines(a,q,ylim = c(0,12),xlim = c(8,10))
abline(v=median(maxnrm$il18prelog))

l<-fitCuminc(time = maxnrm$nrmmn, risk = maxnrm$crn, group =maxnrm$il18prelog>8.876, cens = 0)
o<-l$Tests$stat[2]

o$stat
for( i in (2:588)){
  q[i]<-coxph(Surv(nrmmn,nrmen)~il18prelog>a[i],data=maxnrm)$score
}

plot(a,log(q),ylim = c(0,5),xlim = c(8,10))

abline(v=8.876)


fitCuminc(time = x$nrmmn, risk = x$crn, group = a>8.876, cens = 0)
coxph(Surv(nrmmn,nrmen)~a>8.876,data=x)

summary(coxph(Surv(nrmmn,nrmen)~il18prelog>8.876,data=maxnrm))



h<-coxph(Surv(nrmmn,nrmen)~il18prelog>9,data=maxnrm)


summary(coxph(Surv(nrmmn,nrmen)~maxnrm$il18prelog>maxnrm$il18prelog[20],data=maxnrm))

for( i in (2:588)){
  q[i]<-survfit(time = maxnrm$nrmmn, risk = maxnrm$nrme, group = a>a[i], cens = 0)$Tests[1]$stat[2]
}

b<-survdiff(Surv(nrmmn, nrmen) ~ il18prelog>8.876,data=maxnrm,rho = 0)
summary(b)

fitCuminc(time = maxnrm$nrmmn, risk = maxnrm$crn, group = maxnrm$il18prelog>maxnrm$il18prelog[825], cens = 0)$Tests[1]$stat[2]
p[match(c(i),maxnrm$il18prelog)]<-fitCuminc(time = maxnrm$nrmm, risk = maxnrm$cr, group = groupv, cens = 0)

a<-fitCuminc(time = LUAD$time, risk = LUAD$event, group = LUAD$gender, cens = "alive")

maxnrm$il18prelog[50]
fitCuminc(time = maxnrm$nrmmn, risk = maxnrm$crn, group =(maxnrm$il18prelog>maxnrm$il18prelog[50]), cens = 0)$Tests[1]$stat[2]

a<-etmCIF(Surv(nrmmn,crn>0)~maxnrm$il18prelog>maxnrm$il18prelog[50], etype =crn, data=maxnrm)
maxnrm$crn
summary(a)
  
c<-coxph(Surv(nrmmn,nrmen)~il18prelog>8.876,data=maxnrm)
summary(c)

###########################ADD################################
setwd('D:/work/work table/IL18')
ilhd <- mutate(tbl_df(read.csv('IL18HDadd.CSV')),il18prelog=log2(freeIL18pre),age=Age,cr1=as.factor(cr),freeIL18d28=as.numeric(freeIL18d28),
               Gratwohl2vsrest=Gratwohl2vsrest==2, mismatch3rest1=mismatch3rest1==3, donorsex=donorsex==2,sex=sex==2,myeloid1lymphoid2=myeloid1lymphoid2==2,
               atg=atg==1,ric2rest1=ric2rest1==2,IL18BPaprelog=log2(IL18BPapre), allil18prelog=log2(IL18pre),il18cut=il18prelog>8.876118)
maxnrm<-mutate(ilhd, nrmmn=ifelse(nrmm>60,60,nrmm), nrmen=ifelse(nrmm>60&nrme==1,0,nrme),crn=ifelse(nrmm>60&cr>1,0,cr))%>%filter(!is.na(il18prelog))

coxmax<- coxph(Surv(nrmmn,nrmen)~il18prelog,data=maxnrm)

il18prelog+Age+Karnofsky+BMI30+CVD+Diabetesmellitus+Infection+CRPgreater5+HCTCIkat3
summary(factor(ilhd$Infection))
ggplot(data =filter(ilhd,!is.na(Diagnosis14)), aes(x=factor(Diagnosis14), y=freeIL18pre,na.rm = TRUE)) + geom_boxplot()+theme_bw()  
wilcox.test(freeIL18pre~Karnofsky,data =ilhd)
kruskal.test(freeIL18pre~Diagnosis14,data =ilhd)

ilhd$QAge <- CutQ(ilhd $Age, breaks = quantile(ilhd$Age, seq(1, 0, by = -0.25), na.rm = TRUE), labels = NULL, na.rm = FALSE)
ilhdage<- mutate(ilhd,QAge1=as.numeric(QAge),HCTCI4=plyr::mapvalues(HCTCI,from = c(0:10),to=c(0,1,1,2,2,3,3,3,3,3,3)) )
ggplot(data =filter(ilhdage,!is.na(HCTCI4)), aes(x=factor(HCTCI4), y=freeIL18pre,na.rm = TRUE)) + geom_boxplot()+theme_bw()

jonckheere.test(ilhdage$freeIL18pre, ilhdage$HCTCI4, nperm=5000)
0:10

maxnrm<-mutate(ilhd, nrmmn=ifelse(nrmm>60,60,nrmm), nrmen=ifelse(nrmm>60&nrme==1,0,nrme),crn=ifelse(nrmm>60&cr>1,0,cr))%>%filter(!is.na(il18prelog),Diagnosis14!=4)
coxmax<- coxph(Surv(nrmmn,nrmen)~il18prelog,data=maxnrm)
summary(coxmax)

ilhdossplit <- survSplit(Surv(osm, ose) ~ ., data= ilhd, cut=c(60,125),
                         episode= "tgroup", id="id")
ilhdos1y <- mutate(filter(ilhdossplit,tgroup==1))
ilhdos1ya <- mutate(ilhd, osmn=ifelse(osm>=60,60,osm),osen=ifelse(osm>60,0,ose))%>%filter(Diagnosis14!=4)
hdnrmm1y.cox<- coxph(Surv(osm, ose)~il18prelog,data=ilhdos1ya)
hdnrmm1y.cox<- coxph(Surv(osm, ose)~il18prelog+Age+Karnofsky+BMI30+CRPgreater5+HCTCIkat3,data=ilhdos1y)
summary(hdnrmm1y.cox)

ilhdnrmmsplit <- survSplit(Surv(nrmm, nrme) ~ ., data= ilhd, cut=c(60,100),
                           episode= "tgroup", id="id")
ilhdnrmm1y <-  mutate(filter(ilhdnrmmsplit,tgroup==1))
hdnrmm1y.cox<- coxph(Surv(nrmm, nrme)~il18prelog,data=ilhdnrmm1y)
hdnrmm1y.cox<- coxph(Surv(nrmm, nrme)~il18prelog+Age+Karnofsky+BMI30+CRPgreater5+HCTCIkat3,data=ilhdnrmm1y)
summary(hdnrmm1y.cox)

ilhdnrmmsplit <- survSplit(Surv(nrmm, relapse) ~ ., data= ilhd, cut=c(60,100),
                           episode= "tgroup", id="id")
ilhdnrmm1y <-  mutate(filter(ilhdnrmmsplit,tgroup==1))
hdnrmm1y.cox<- coxph(Surv(nrmm, relapse)~il18prelog+Age+Karnofsky+BMI30+CRPgreater5+HCTCIkat3,data=ilhdnrmm1y)
summary(hdnrmm1y.cox)



