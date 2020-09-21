#########################################################################################
### FUNKTIONEN #########################################################################

maxgray = function(time,status, cov, alpha = NA, minprob = 0.1, maxprob = 0.9, event = 1, cens = 0, compevent = 2){
  
  # Datenvorbereitung

  data = data.frame(time,status,cov)
 data = data[complete.cases(data[,3]),]
  data = data[order(cov),]  
  data$status =ifelse(data$status == cens, 0, ifelse(data$status ==event, 1,2))
  
  # Selektion der Cutpoints
  if(minprob < 0 || maxprob > 1) stop("minprob and/or maxprob not in intervall [0,1]")
  if( minprob > maxprob) stop("minprob < maxprob not given")
  m = quantile(data$cov, c(minprob,maxprob))
  work = data
  get.cuts = which(work$cov > m[1] & (work$cov < m[2]))
  if(length(get.cuts) < 1) stop("no data between minprob and maxprob")
  cuts = unique(work$cov[get.cuts])
  
  # Berechnung der selektierten Teststatistiken
  tests = sapply(cuts, graystatm, data=work)
  
  #maximal selektierte Teststatistik
  maxtest = max(tests)
  get.max = which( maxtest ==tests)
  
  #optimaler Cutpoint
  cut = unique(cuts)[min(get.max)]
  
  #P-Wert der maximale selektierten Teststatistik
  epsilon = 0.0001
  p1 = 0
  i = 1
  repeat{
    p2 = p1+((-1)^(i+1))*exp((-2)*(i^2)*(maxtest^2))
    if(abs(p2-p1) < epsilon/2) {break}
    p1 = p2
    i= i+1
  }
  pval = p2*2
  
  # Quantil
  if(!is.null(alpha)) {quant = sqrt(-(0.5*log(alpha/2)))}
  if (is.null(alpha)) {quant <- NA}
  
  #Ausgabe
  #pars = as.list(match.call())
  pars = as.character(match.call())
  pars = unlist(strsplit(pars[4], "$", fixed = TRUE))
  sumevents = nrow(subset(data, status ==1))
  sumcens = nrow(subset(data, status ==0))
  sumcr = nrow(subset(data, status ==2))
  out.list = list(maxtest = maxtest, estimate = cut, stats = tests, cuts = cuts, pvalue = pval,
                  nevents = sumevents, ncens = sumcens, ncompev = sumcr, quant = quant, factor = pars[2]) # factor= as.character(pars$cov)
  class(out.list) = "maxgray"
  out.list
}


graystatm = function(cuts, data){
  data$group = NA
  for( i in 1:nrow(data)){
    if(data$cov[i] < cuts) {data$group[i] = 0} 
    if(data$cov[i] >= cuts) {data$group[i] = 1} }
  
  tests1 = NA
  tests1 = graystat(data$time, data$status, data$group)
  return(tests1)
}


graystat = function(time, status, group){
  
  if ("reshape" %in% rownames(installed.packages())==F) stop( "Require package reshape")
  data = data.frame(time,status,group)
  
  #Datenvorbereitung  
  data$group = data$group +1
  temp = data
  ng = 2
  dat.long = data[,c("time", "status", "group")]
  
  # Variablen für Ereignizeiten nach Status und Gruppe
  dat.long$freq = 1
  dat.long$idnew = c(rep(0, nrow(dat.long)))
  dat.long = dat.long[order(dat.long$time, dat.long$status, dat.long$group, dat.long$freq),]  
  
  # Ereignis d_status_group
  for (i in 1:nrow(dat.long)){
    dat.long[i,5] = paste("d_",dat.long[i,2],"_",dat.long[i,3], sep="")
  }
  dat.long$temp = c(1:nrow(dat.long))
  require(reshape)
  dat.long = melt(dat.long, c("time","idnew", "temp"), "freq")
  dat.long = cast(dat.long, temp+time~idnew)
  dat.long[is.na(dat.long)] = 0
  
  if(("d_0_1" %in% colnames(dat.long)) == FALSE) { dat.long$d_0_1 = 0}
  if(("d_0_2" %in% colnames(dat.long)) == FALSE) { dat.long$d_0_2 = 0}
  if(("d_1_1" %in% colnames(dat.long)) == FALSE) { dat.long$d_1_1 = 0}
  if(("d_1_2" %in% colnames(dat.long)) == FALSE) { dat.long$d_1_2 = 0}
  if(("d_2_1" %in% colnames(dat.long)) == FALSE) { dat.long$d_2_1 = 0}
  if(("d_2_2" %in% colnames(dat.long)) == FALSE) { dat.long$d_2_2 = 0}
  
  #Datensatz mit Ereignisvariablen
  dat.long = dat.long[,c( "time", "d_0_1", "d_0_2", "d_1_1", "d_1_2", "d_2_1", "d_2_2")]
  
  # Hinzufügen erste Zeile zum Datensatz
  firstrow = as.data.frame(t(c(rep(0, ncol(dat.long)))))
  colnames(firstrow) = colnames(dat.long)
  dat.long = rbind(firstrow, dat.long)
  row.names(dat.long) = 0:(nrow(dat.long)-1)
  
  # Initialisierung der Variablen
  dat.long$nd1 = dat.long$d_1_1 + dat.long$d_1_2 # interessierende Ereignisse 
  dat.long$nd2 = dat.long$d_2_1 + dat.long$d_2_2 # konkurrierende Ereignisse
  dat.long$td1 = dat.long$d_1_1 + dat.long$d_2_1 # Ereignisse Gruppe 1
  dat.long$td2 = dat.long$d_1_2 + dat.long$d_2_2 # Ereignisse Gruppe 2
  dat.long$rs_1 = nrow(subset(temp, group ==1)) # Anzahl unter Risiko Gruppe 1
  dat.long$rs_2 = nrow(subset(temp, group ==2)) # Anzahl unter Risiko Gruppe 2
  dat.long$rsh1 = 1 #Hilfsvariable
  dat.long$rsh2 = 1 # Hilfsvariable
  dat.long$skm_1 = 1 #Überlebensfunktion Gruppe 1
  dat.long$skm_2 = 1 #Überlebensfunktion Gruppe 2
  dat.long$f1_1 = 0 # Verteilungsfunktion Gruppe 1
  dat.long$f1_2 = 0 # Verteilungsfunktion Gruppe 2
  dat.long$weight_1 = 0 # Gewicht Gruppe 1
  dat.long$weight_2 = 0 # Gewicht Gruppe 2
  dat.long$rw_1 = 0 # adjustiertes Risiko Gruppe 1 (Gewicht*Risiko)
  dat.long$rw_2 = 0 # adjustiertes Risiko Gruppe 2 (Gewicht*Risiko)
  
  get.rs_ = which(grepl("rs_", colnames(dat.long)) ==TRUE)
  get.rw = which(grepl("rw", colnames(dat.long)) ==TRUE)
  get.nd = which(grepl("nd", colnames(dat.long)) ==TRUE)
  get.d_0 = which(grepl("d_0_", colnames(dat.long)) == TRUE)
  get.d_1 = which(grepl("d_1_", colnames(dat.long)) == TRUE)
  get.f1_ = which(grepl("f1_", colnames(dat.long)) == TRUE)
  get.skm_ = which(grepl("skm_", colnames(dat.long)) == TRUE)
  get.td = which(grepl("td", colnames(dat.long)) == TRUE)
  get.tr = which(grepl("tr", colnames(dat.long)) == TRUE)
  get.weight_ = which(grepl("weight_", colnames(dat.long)) == TRUE)
  get.rsh = which(grepl("rsh", colnames(dat.long)) ==TRUE)
  
  for ( j in 1:ng){
    for ( i in 2:nrow(dat.long)) {
      dat.long[i,get.rs_[j]] = dat.long[i-1, get.rs_[j]] - dat.long[i, get.td[j]] - dat.long[i, get.d_0[j]]
    }}
  
  for ( j in 1:ng){
    for ( i in 2:nrow(dat.long)) {
      dat.long[i,get.rsh[j]] = dat.long[i-1,get.rs_[j]]
    }}
  
  #Berechnung weight = (1-F)/S nach Grays-Test 2.2.2
  for ( j in 1:ng){
    for ( i in 2:nrow(dat.long)) {
      
      if(dat.long[i, "nd1"] || dat.long[i, "nd2"]  > 0) {
        if(dat.long[i,get.rsh[j]] > 0){
          dat.long[i, get.skm_[j]] = dat.long[i-1, get.skm_[j]] * ( dat.long[i-1, get.rs_[j]] - dat.long[i, get.td[j]])/ dat.long[i-1, get.rs_[j]]
          dat.long[i,get.f1_[j]] =  dat.long[i-1,get.f1_[j]] + (dat.long[i-1,get.skm_[j]]* dat.long[i,get.d_1[j]]) / dat.long[i-1,get.rs_[j]]
          dat.long[i, get.weight_[j]] = (1 - dat.long[i-1, get.f1_[j]]) / dat.long[i-1, get.skm_[j]]
          dat.long[i, get.rw[j]] = dat.long[i-1, get.rs_[j]] * dat.long[i, get.weight_[j]] } 
        else{dat.long[i, get.skm_[j]]= 0
        dat.long[i, get.f1_[j]] = dat.long[i-1, get.f1_[j]]}}
      else{dat.long[i, get.skm_[j]]= dat.long[i-1, get.skm_[j]]
      dat.long[i, get.f1_[j]]= dat.long[i-1, get.f1_[j]]}
    } }
  
  #Score
  aa = dat.long[ -1,]
  aa$rs_w = aa$rw_1 + aa$rw_2 # adjustiertes Risiko total
  aa$eta = NA
  aa$group_00 = aa$d_0_1 + aa$d_1_1 + aa$d_2_1 #Beobachtungen Gruppe 1
  aa$group_11 = aa$d_0_2 + aa$d_1_2 + aa$d_2_2 #Beobachtungen Gruppe 2
  aa$weight0_R = 0
  aa$weight1_R = 0
  
  for (i in 1:nrow(aa)){
    if(aa[i, "nd1"] >= 1){ aa[i, "eta"] = 1} else {aa[i, "eta"] = 0}
    if( aa$rs_w[i] > 0) {
      aa$weight0_R[i] = aa[i, "weight_1"] / aa[i, "rs_w"] * aa[i, "nd1"]
      aa$weight1_R[i] = aa[i, "weight_2"] / aa[i, "rs_w"] * aa[i, "nd1"]
    }
  }
  
  bb = aa[, c("eta", "group_00", "group_11", "weight0_R", "weight1_R")]
  bb$sum_weight0 = cumsum(bb$weight0_R)
  bb$sum_weight1 = cumsum(bb$weight1_R)
  bb$sum_w0_r = bb$sum_weight0 * bb$group_00
  bb$sum_w1_r = bb$sum_weight1 * bb$group_11
  bb$a = bb$sum_w0_r + bb$sum_w1_r
  bb$s = bb$eta - bb$a
  bb$s2 = bb$s ^ 2
  G = abs(sum(bb$s * bb$group_11))
  sum_s2 = sum(bb$s2)
  sqrt_sum_s2=sqrt(sum_s2)
  #Teststatistik
  Q = G/sqrt_sum_s2
  
  return(Q) 
}



print.maxgray = function(x, digits = getOption("digits"), ...) {
  cat("\n Maximally selected Gray Statistic \n for competing risks \n")
  cat("\n")
  cat("number of events of interest = ", x$nevents, "\n")
  cat("number of competing events = ", x$ncompev, "\n")
  cat("number of censorings = ", x$ncens, "\n")
  cat("M = ", x$maxtest, " , p-value = ", x$pvalue, "\n")
  cat("estimated cutpoint = ", x$estimate, "\n")
  
  #cat( x$cuts, "\n", x$stats )
}


plot.maxgray = function(x,xlab = NULL, ylab = NULL, ...){
  
  if(is.null(ylab)) { ylab = "Standardized Gray's-Statistic"}
  if(is.null(xlab)) {xlab = x$factor}
  if(!is.na(x$quant)){
    ylim <- c(min(x$quant, min(x$stats)), max(x$quant, max(x$stats)))
    ylim <- c(ylim[1]*0.95, ylim[2]*1.05)
    xlength = range(x$cuts)
    plot(x$cuts, x$stats, type = "b", xlab = xlab, ylab = ylab, ylim = ylim,...)
    lines(c(x$estimate, x$estimate), c(0,x$maxtest), lty = 3)
    lines(xlength, c(x$quant, x$quant), col=2)
  }
  if(is.na(x$quant)){
    xlength = range(x$cuts)
    plot(x$cuts, x$stats, type = "b", xlab = xlab, ylab = ylab,...)
    lines(c(x$estimate, x$estimate), c(0,x$maxtest), lty = 3)
  }
}


######################################################################################################################
##### ANWENDUNG  #################################################################################################

load("C:/Users/Hannah Schmidt/Documents/Bachelorarbeit/Daten/VitD.RData")
data = dat.train[, c("vitD", "TimetoRelapse", "CompStatusnum")]
data$CompStatusnum = as.factor(data$CompStatusnum)

###########################################################################################################################################
### Competing Risk Cutpoint
fm1 = maxgray(data$TimetoRelapse, data$CompStatusnum, data$vitD, alpha = 0.05, event =3, cens = 1, compevent = 2)
print(fm1)
plot(fm1)
index<-function(x,data){
  data[[x]]
}
index(nrmm,vodh)

