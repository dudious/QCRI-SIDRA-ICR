### R code from vignette source 'vignettes/genefu/inst/doc/genefu.Rnw'

###################################################
### code chunk number 1: loadPackages
###################################################
library(genefu)
library(xtable)
library(rmeta)


###################################################
### code chunk number 2: installAndLoadMAINZ (eval = FALSE)
###################################################
## source("http://www.bioconductor.org/biocLite.R")
## biocLite("breastCancerMAINZ")
## biocLite("breastCancerTRANSBIG")
## biocLite("breastCancerUPP")
## biocLite("breastCancerUNT")
## biocLite("breastCancerNKI")


###################################################
### code chunk number 3: installAndLoadAllPackages
###################################################
library(breastCancerMAINZ)
library(breastCancerTRANSBIG)
library(breastCancerUPP)
library(breastCancerUNT)
library(breastCancerNKI)


###################################################
### code chunk number 4: findDuplicatedPatients
###################################################
library(Biobase)
data(breastCancerData)
cinfo <- colnames(pData(mainz7g))
data.all <- c("transbig7g"=transbig7g, "unt7g"=unt7g, "upp7g"=upp7g, "mainz7g"=mainz7g, "nki7g"=nki7g)

idtoremove.all <- NULL
duplres <- NULL

## UNT vs UPP vs TRANSBIG
demo.all <- rbind(pData(transbig7g), pData(unt7g), pData(upp7g))
dn2 <- c("TRANSBIG", "UNT", "UPP")

## Karolinska
## VDXKIU, KIU, UPPU
ds2 <- c("VDXKIU", "KIU", "UPPU")
demot <- demo.all[complete.cases(demo.all[ , c("series")]) & is.element(demo.all[ , "series"], ds2), ]

duplid <- sort(unique(demot[duplicated(demot[ , "id"]), "id"]))
duplrest <- NULL
for(i in 1:length(duplid)) {
  tt <- NULL
  for(k in 1:length(dn2)) {
    myx <- sort(row.names(demot)[complete.cases(demot[ , c("id", "dataset")]) & demot[ , "id"] == duplid[i] & demot[ , "dataset"] == dn2[k]])
    if(length(myx) > 0) { tt <- c(tt, myx) }
  }
  duplrest <- c(duplrest, list(tt))
}
names(duplrest) <- duplid
duplres <- c(duplres, duplrest)

## Oxford
## VVDXOXFU, OXFU
ds2 <- c("VDXOXFU", "OXFU")
demot <- demo.all[complete.cases(demo.all[ , c("series")]) & is.element(demo.all[ , "series"], ds2), ]

duplid <- sort(unique(demot[duplicated(demot[ , "id"]), "id"]))
duplrest <- NULL
for(i in 1:length(duplid)) {
  tt <- NULL
  for(k in 1:length(dn2)) {
    myx <- sort(row.names(demot)[complete.cases(demot[ , c("id", "dataset")]) & demot[ , "id"] == duplid[i] & demot[ , "dataset"] == dn2[k]])
    if(length(myx) > 0) { tt <- c(tt, myx) }
  }
  duplrest <- c(duplrest, list(tt))
}
names(duplrest) <- duplid
duplres <- c(duplres, duplrest)

## duplicated patients
duPL <- sort(unlist(lapply(duplres, function(x) { return(x[-1]) } )))


###################################################
### code chunk number 5: computeRiskScore
###################################################
dn <- c("transbig", "unt", "upp", "mainz", "nki")
dn.platform <- c("affy", "affy", "affy", "affy", "agilent")

res <- ddemo.all <- ddemo.coln <- NULL
for(i in 1:length(dn)) {

  ## load dataset
  dd <- get(data(list=dn[i]))

  ddata <- t(exprs(dd))
  ddemo <- phenoData(dd)@data
  dannot <- featureData(dd)@data
  ddemo.all <- c(ddemo.all, list(ddemo))
  if(is.null(ddemo.coln)) { ddemo.coln <- colnames(ddemo) } else { ddemo.coln <- intersect(ddemo.coln, colnames(ddemo)) }
  rest <- NULL

  ## NPI
  ss <- ddemo[ , "size"]
  gg <- ddemo[ , "grade"]
  nn <- rep(NA, nrow(ddemo))
  nn[complete.cases(ddemo[ , "node"]) & ddemo[ , "node"] == 0] <- 1 
  nn[complete.cases(ddemo[ , "node"]) & ddemo[ , "node"] == 1] <- 3
  names(ss) <- names(gg) <- names(nn) <- rownames(ddemo)
  rest <- cbind(rest, "NPI"=npi(size=ss, grade=gg, node=nn, na.rm=TRUE)$score)

  ## AURKA
  ## if affy platform consider the probe published in Desmedt et al., CCR, 2008
  if(dn.platform[i] == "affy") { domap <- FALSE } else { domap <- TRUE }
  modt <- scmgene.robust$mod$AURKA
  ## if agilent platform consider the probe published in Desmedt et al., CCR, 2008
  if(dn.platform[i] == "agilent") {
    domap <- FALSE
    modt[ , "probe"] <- "NM_003600"
  }
  rest <- cbind(rest, "AURKA"=sig.score(x=modt, data=ddata, annot=dannot, do.mapping=domap)$score)

  ## GGI
  if(dn.platform[i] == "affy") { domap <- FALSE } else { domap <- TRUE }
  rest <- cbind(rest, "GGI"=ggi(data=ddata, annot=dannot, do.mapping=domap)$score)
  ## GENIUS
  if(dn.platform[i] == "affy") { domap <- FALSE } else { domap <- TRUE }
  rest <- cbind(rest, "GENIUS"=genius(data=ddata, annot=dannot, do.mapping=domap)$score)
  res <- rbind(res, rest)
}
names(ddemo.all) <- dn


###################################################
### code chunk number 6: simplifyAndRemoveDuplicatePatients
###################################################
ddemot <- NULL
for(i in 1:length(ddemo.all)) {
  ddemot <- rbind(ddemot, ddemo.all[[i]][ , ddemo.coln, drop=FALSE])
}
res[complete.cases(ddemot[ ,"dataset"]) & ddemot[ ,"dataset"] == "VDX", "GENIUS"] <- NA

## select only untreated node-negative patients with all risk predictions
myx <- complete.cases(res, ddemot[ , c("node", "treatment")]) & ddemot[ , "treatment"] == 0 & ddemot[ , "node"] == 0 & !is.element(rownames(ddemot), duPL)

res <- res[myx, , drop=FALSE]
ddemot <- ddemot[myx, , drop=FALSE]


###################################################
### code chunk number 7: cindexComputation
###################################################
cc.res <- complete.cases(res)
datasetList <- c("MAINZ","TRANSBIG","UPP","UNT","NKI")
riskPList <- c("NPI", "AURKA", "GGI", "GENIUS")
setT <- setE <- NULL
resMatrix <- as.list(NULL)

for(i in datasetList){
  dataset.only <- ddemot[,"dataset"] == i
  patientsAll <- cc.res & dataset.only
	
  ## set type of available survival data
  if(i == "UPP") {
    setT <- "t.rfs"
    setE <- "e.rfs"
  } else {
    setT <- "t.dmfs"
    setE <- "e.dmfs"
  }

  ## cindex computation
  cindexN <- t(apply(X=t(res[patientsAll,"NPI"]), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=ddemot[patientsAll,setT], z=ddemot[patientsAll, setE]))
	
  cindexA <- t(apply(X=t(res[patientsAll,"AURKA"]), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=ddemot[patientsAll,setT], z=ddemot[patientsAll, setE]))
	
  cindexM <- t(apply(X=t(res[patientsAll,"GGI"]), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=ddemot[patientsAll, setT], z=ddemot[patientsAll, setE]))
	
  cindexG <- t(apply(X=t(res[patientsAll,"GENIUS"]), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=ddemot[patientsAll, setT], z=ddemot[patientsAll, setE]))

	
  resMatrix[["NPI"]] <- rbind(resMatrix[["NPI"]], cindexN)
  resMatrix[["AURKA"]] <- rbind(resMatrix[["AURKA"]], cindexA)
  resMatrix[["GGI"]] <- rbind(resMatrix[["GGI"]], cindexM)
  resMatrix[["GENIUS"]] <- rbind(resMatrix[["GENIUS"]], cindexG)
}


###################################################
### code chunk number 8: combineEstimations
###################################################
for(i in names(resMatrix)){
  ceData <- combine.est(x=resMatrix[[i]][,"cindex"], x.se=resMatrix[[i]][,"cindex.se"], hetero=TRUE)
  cLower <- ceData$estimate + qnorm(0.025, lower.tail=TRUE) * ceData$se
  cUpper <- ceData$estimate + qnorm(0.025, lower.tail=FALSE) * ceData$se
	
  cindexO <- cbind("cindex"=ceData$estimate, "cindex.se"=ceData$se, "lower"=cLower, "upper"=cUpper)
  resMatrix[[i]] <- rbind(resMatrix[[i]], cindexO)
  rownames(resMatrix[[i]]) <- c(datasetList, "Overall")
}


###################################################
### code chunk number 9: computePValues
###################################################
pv <- sapply(resMatrix, function(x) { return(x["Overall", c("cindex","cindex.se")]) })
pv <- apply(pv, 2, function(x) { return(pnorm((x[1] - 0.5) / x[2], lower.tail=x[1] < 0.5)) })
printPV <- matrix(pv,ncol=4)
rownames(printPV) <- "P-value"
colnames(printPV) <- names(pv)


###################################################
### code chunk number 10: printPvalue
###################################################
xtable(printPV, digits=c(0, rep(-1,ncol(printPV))))


###################################################
### code chunk number 11: forestplotNPI
###################################################
par(mfrow=c(2,2))
datasetListF <- c("MAINZ","TRANSBIG","UPP","UNT","NKI", NA, "Overall")
myspace <- "   "

## NPI Forestplot
tt <- rbind(resMatrix[["NPI"]][1:5,],
          "NA"=NA,
          "Overall"=resMatrix[["NPI"]][6,])

tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset", datasetListF))

r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)

metaplot.surv(mn=r.mean, lower=r.lower, upper=r.upper, labels=labeltext, xlim=c(0.4,0.9), boxsize=0.5, zero=0.5, col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"), main="NPI Concordance Index")
#@
#
#<<forestplotAURKA,fig=TRUE,echo=FALSE>>=
## AURKA Forestplot
tt <- rbind(resMatrix[["AURKA"]][1:5,],
          "NA"=NA,
          "Overall"=resMatrix[["AURKA"]][6,])

tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset", datasetListF))

r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)

metaplot.surv(mn=r.mean, lower=r.lower, upper=r.upper, labels=labeltext, xlim=c(0.4,0.9), boxsize=0.5, zero=0.5, col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"), main="AURKA Concordance Index")
#@
#
#<<forestplotGGI,fig=TRUE,echo=FALSE>>=
## GGI Forestplot
tt <- rbind(resMatrix[["GGI"]][1:5,],
          "NA"=NA,
          "Overall"=resMatrix[["GGI"]][6,])

tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset", datasetListF))

r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)

metaplot.surv(mn=r.mean, lower=r.lower, upper=r.upper, labels=labeltext, xlim=c(0.4,0.9), boxsize=0.5, zero=0.5, col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"), main="GGI Concordance Index")
#@
#
#<<forestplotGENIUS,fig=TRUE,echo=FALSE>>=
## GENIUS Forestplot
tt <- rbind(resMatrix[["GENIUS"]][1:5,],
          "NA"=NA,
          "Overall"=resMatrix[["GENIUS"]][6,])

tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset", datasetListF))

r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)

metaplot.surv(mn=r.mean, lower=r.lower, upper=r.upper, labels=labeltext, xlim=c(0.4,0.9), boxsize=0.5, zero=0.5, col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"), main="GENIUS Concordance Index")


###################################################
### code chunk number 12: forestplotOverall
###################################################
## Overall Forestplot
mybigspace <- "       "
tt <- rbind("OverallN"=resMatrix[["NPI"]][6,],
          "OverallA"=resMatrix[["AURKA"]][6,],
          "OverallM"=resMatrix[["GGI"]][6,],
          "OverallG"=resMatrix[["GENIUS"]][6,])

tt <- as.data.frame(tt)
labeltext <- cbind(c("Risk Prediction","NPI","AURKA","GGI","GENIUS"))

r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)

metaplot.surv(mn=r.mean, lower=r.lower, upper=r.upper, labels=labeltext, xlim=c(0.45,0.75), boxsize=0.5, zero=0.5, col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"), main="Overall Concordance Index")


###################################################
### code chunk number 13: computeCindexWithPvalue
###################################################
cc.res <- complete.cases(res)
datasetList <- c("MAINZ","TRANSBIG","UPP","UNT","NKI")
riskPList <- c("NPI", "AURKA", "GGI", "GENIUS")
setT <- setE <- NULL
resMatrixFull <- as.list(NULL)

for(i in datasetList){
  dataset.only <- ddemot[,"dataset"] == i
  patientsAll <- cc.res & dataset.only
	
  ## set type of available survival data
  if(i == "UPP") {
    setT <- "t.rfs"
    setE <- "e.rfs"
  } else {
    setT <- "t.dmfs"
    setE <- "e.dmfs"
  }
	
  ## cindex and p-value computation
  cindexN <- t(apply(X=t(res[patientsAll,"NPI"]), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(tt); },
    y=ddemot[patientsAll,setT], z=ddemot[patientsAll, setE]))
	
  cindexA <- t(apply(X=t(res[patientsAll,"AURKA"]), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(tt); },
    y=ddemot[patientsAll,setT], z=ddemot[patientsAll, setE]))
	
  cindexM <- t(apply(X=t(res[patientsAll,"GGI"]), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(tt); },
    y=ddemot[patientsAll, setT], z=ddemot[patientsAll, setE]))
	
  cindexG <- t(apply(X=t(res[patientsAll,"GENIUS"]), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(tt); },
    y=ddemot[patientsAll, setT], z=ddemot[patientsAll, setE]))
	
  resMatrixFull[["NPI"]] <- rbind(resMatrixFull[["NPI"]], cindexN)
  resMatrixFull[["AURKA"]] <- rbind(resMatrixFull[["AURKA"]], cindexA)
  resMatrixFull[["GGI"]] <- rbind(resMatrixFull[["GGI"]], cindexM)
  resMatrixFull[["GENIUS"]] <- rbind(resMatrixFull[["GENIUS"]], cindexG)
}

for(i in names(resMatrixFull)){
  rownames(resMatrixFull[[i]]) <- datasetList
}

ccmData <- tt <- rr <- NULL
for(i in 1:length(resMatrixFull)){
  tt <- NULL
  for(j in 1:length(resMatrixFull)){
    if(i != j) { rr <- cindex.comp.meta(list.cindex1=resMatrixFull[[i]], list.cindex2=resMatrixFull[[j]], hetero=TRUE)$p.value } else { rr <- 1 }
    tt <- cbind(tt, rr)
  }
  ccmData <- rbind(ccmData, tt)
}
ccmData <- as.data.frame(ccmData)
colnames(ccmData) <- riskPList
rownames(ccmData) <- riskPList


###################################################
### code chunk number 14: printCCM
###################################################
xtable(ccmData, digits=c(0, rep(-1,ncol(ccmData))))


###################################################
### code chunk number 15: computeCCMPval
###################################################
ccmDataPval <- matrix(p.adjust(data.matrix(ccmData), method="holm"),ncol=4,dimnames=list(rownames(ccmData),colnames(ccmData)))


###################################################
### code chunk number 16: printCCMPval
###################################################
xtable(ccmDataPval, digits=c(0, rep(-1,ncol(ccmDataPval))))


###################################################
### code chunk number 17: sessionInfo
###################################################
toLatex(sessionInfo())


