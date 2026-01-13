##  Set a working directory
setwd("/home/leej40/Documents/extlinear/data")
load("cumbriaLanc.RData")
#should load four R objects: prcp, stations, dateVec, dat
#prcp is dat without first column of dates to make indexing more simple
prcp=as.data.frame(prcp)
stations=as.data.frame(stations)
dat=as.data.frame(dat)

##  EDA
numDays <- dim(prcp)[1]  #[1] 56274
numStations <- dim(stations)[1]  #[1] 608
pctMissing <- apply(is.na(prcp), 2, sum)/numDays
# retaining data only after 1960
date1960 <- which(dateVec == "1960-01-01")  #32506
prcp <- prcp[date1960:numDays,]
dat <- dat[date1960:numDays,]
dateVec <- dateVec[date1960:numDays]
numDays <- dim(prcp)[1]
pctMissing <- apply(is.na(prcp), 2, sum)/numDays

## retaining stations which have more missing < .1
keep <- pctMissing < .1
sum(keep)
prcp <- prcp[,keep]
dat <- dat[, c(1,which(keep)+1)]
stations <- stations[keep,]
numStations <- dim(stations)[1] #30
pctMissing <- pctMissing[keep]
prcp=na.omit(prcp)
dim(prcp)
dat=na.omit(dat)
dim(dat)

##  find max of each station
maxAtStat <- apply(prcp, 2, max, na.rm = T)
range(maxAtStat)
##  look at day which is has absolute maximum as sanity check
statWithMax <- which(maxAtStat == max(maxAtStat)) #66
dayOfMax <- which(prcp[, statWithMax] == max(maxAtStat)) #3375
maxDayPrcp <- prcp[dayOfMax,]
#keep <- !is.na(maxDayPrcp)
keep <- as.vector(!is.na(maxDayPrcp))
library(fields)
#ai <- as.image(as.vector(maxDayPrcp[keep]), x = stations[keep, c("lon", "lat")])
# Convert the 1-row tibble to a numeric vector
vals <- as.numeric(maxDayPrcp[keep])
# Convert the coordinate tibble to a matrix
coords <- as.matrix(stations[keep, c("lon", "lat")])
# Run as.image with the correct data types
ai <- as.image(vals, x = coords)
dev.new()
image.plot(ai)

#look at xi estimates
xiEst <- numeric(numStations)
library(lmomco)
for(i in seq(1, numStations)){
  prcpi <- prcp[!is.na(prcp[,i]),i]
  q.99 <- quantile(prcpi, .99)
  prcpiExc <- prcpi[prcpi > q.99]
  lmomEsts <- lmoms(prcpiExc)
  xiEst[i] <- -pargpa(lmomEsts)$para[3]
}  
mean(xiEst) #0.064
mean(xiEst>0)  #0.73 have positive estimates
dev.new()
plot(xiEst)  
abline(h = 0, lty = 3)
#some are pretty negative!

##  Divide the data set into two sets randomly
##  Rank transform to Pareto-2 marginals
set.seed(12345)
n_obs=dim(prcp)[1]
R_index=sample(x = seq(1:n_obs),size = n_obs,replace = F)
prcp=prcp[R_index,]
#prcp=as.data.frame(prcp)
trun=ceiling(n_obs*(2/3))
prcp_tr=prcp[1:trun,]
prcp_te=prcp[-(1:trun),]

U_tr <- apply(prcp_tr, 2, function(i) rank(i)/(dim(prcp_tr)[1]+1))
U_te <- apply(prcp_te, 2, function(i) rank(i)/(dim(prcp_te)[1]+1))
##  Simulate from a shifted Pareto distribution
shift=0.9352074
fPrcp_tr <- apply(U_tr, 2, function(x) (1-x)^(-1/2)-shift)
fPrcp_te <- apply(U_te, 2, function(x) (1-x)^(-1/2)-shift)

##  Pairwise TPDM estimation
TPDM <- matrix(nrow = numStations, ncol = numStations)
for(i in seq(1, numStations-1)){
  for(j in seq(i+1, numStations)){
    fPrcp_ij <- fPrcp_tr[, c(i, j)]
    rad <- sqrt(apply(fPrcp_ij^2, 1, sum))
    ang <- fPrcp_ij/rad
    q.99 <- quantile(rad, 0.95)
    radKeep <- rad > q.99
    angKeep <- ang[radKeep,]
    tpdEst <- 2*mean(apply(angKeep, 1, prod))
    TPDM[i,j] <- TPDM[j,i] <- tpdEst
  }
}
diag(TPDM) <- 1
sum(is.na(TPDM))  #0

#makes a map showing the strength of dependence (by TPDM)
#with location 2, chosen as it has lowest pctMissing
TPDM_store <- TPDM  #cycle thru TPDM's
loi <- which(pctMissing == min(pctMissing))  #4 if you don't want to mess with reordering
stations_mtx <- as.matrix(stations[, c("lon", "lat")])
ai <- as.image(TPDM_store[loi,], x = stations_mtx)
pdf("/home/leej40/Documents/extlinear/TPDM_final.pdf",6,6)
dev.new()
par(mar = c(5, 5, 4, 1))
image.plot(ai,cex.main=1.5, cex.lab=1.5, cex.axis=1.5, cex=1.5,axis.args = list(cex.axis = 1.2),xlab="lon",ylab="lat")
points(stations[loi, "lon"]+0.01, stations[loi, "lat"],pch=8, cex=1,lwd=1)
dev.off()

# what would prediction weights be if all stations were available?
TPDM_store <- TPDM 
idx <- 1:numStations
idx <- c(idx[-loi], idx[loi])
TPDM_reord <- TPDM_store[idx, idx]
Sig_11 <- TPDM_reord[-numStations, -numStations]
Sig_12 <- TPDM_reord[-numStations, numStations]
b <- solve(Sig_11, Sig_12)
stations_mtx <- as.matrix(stations[-loi, c("lon", "lat")])
ai <- as.image(b, x = stations_mtx)
pdf("/home/leej40/Documents/extlinear/weights_final.pdf",6,6)
dev.new()
par(mar = c(5, 5, 4, 1))
image.plot(ai, cex.main=1.5, cex.lab=1.5, cex.axis=1.5, cex=1.3, axis.args = list(cex.axis = 1.2),xlab="lon",ylab="lat")
points(stations[loi, "lon"]+0.01, stations[loi, "lat"],pch=8, cex=1.1,lwd=1)
dev.off()

##  Prediction procedure
setwd("/home/leej40/Documents/extlinear/code")
source("TransformedOperations.R")
library(MASS) #ginv
source("CPfactor.R")
library(plotrix)  # 'draw.circle'
source("jointRegion.R")
source("crossValidate.R")
source("estimateParams.R")
library(VGAM)
library(ks)
source("KDE_w.R")
source("condDensity.R")
source("coverageRate.R")

#iMtx <- diag(numStations)
numPred <- dim(fPrcp_te)[1]
TPDM_store <- TPDM  #cycle through TPDM's
predVec <- numeric(numPred)
obsVec <- numeric(numPred)

for(i in seq(1:numPred)){
  
  #doiKeep <- !is.na(fPrcp_te[i,])
  #doi_fPrcp <- fPrcp_te[i,doiKeep]
  #doiStations <- stations[doiKeep,]
  #T1 <- iMtx[doiKeep,]
  #TPDM_doi <- T1 %*% TPDM_store %*% t(T1)
  #reorder to put prediction location in last row
  #idx_loi <- which(doiStations$station_id == stations[loi, ]$station_id)
  #P_reduced <- dim(TPDM_doi)[1]
  idx_loi <- which(stations$station_id == stations[loi, ]$station_id)
  P_reduced <- dim(TPDM)[1]
  idx <- seq(1, P_reduced)
  idx <- c(idx[-idx_loi], idx_loi)
  #TPDM_reord <- TPDM_doi[idx, idx]
  TPDM_reord <- TPDM[idx, idx]
  Sig_11 <- TPDM_reord[-P_reduced, -P_reduced]
  Sig_12 <- TPDM_reord[-P_reduced, P_reduced]
  b <- solve(Sig_11, Sig_12)
  #predVec[i] <- t(b) %*% fPrcp_te[-idx_loi]
  predVec <- as.vector(Amul(b, t(fPrcp_te[,-idx_loi])))
  obsVec <- fPrcp_te[,idx_loi]
}
summary(predVec) # Xhat
summary(obsVec) # observations in a test set


idx <- seq(1, dim(TPDM_store)[1])
idx <- c(idx[-loi], loi)
TPDM_reord <- TPDM_store[idx, idx]

d=dim(TPDM_reord)[1]
Sig11=TPDM_reord[1:(d-1),1:(d-1)]
Sig12=TPDM_reord[1:(d-1),d]
Sig21=t(Sig12)
Sig22=TPDM_reord[d,d]
##  Prediction inner product matrix
Sig_P_11=Sig21%*%solve(Sig11)%*%Sig12
TPDM_Phat=matrix(c(Sig_P_11,Sig_P_11,Sig_P_11,Sig22),nrow=2)

####  Decompose the TPDM_Pred or the TPDM_Pred_hat
CP_prcp1=CPfactor(Mtx = TPDM_Phat,q_star = 10,ite_pereach = 5000,ite_cp = 100)
save(CP_prcp1,file="CP_prcp_final.RData")
load(file="CP_prcp_final.RData")

####  Check the coverage rate
####  Create the 95% joint polar region from a CP-factorization
##  Find the 95% joint polar region (Figure 3 (left))
jointOut=jointRegion(Xhat = predVec, Xf = obsVec,
                     Angular = CP_prcp1$angular, Pmass = CP_prcp1$pmass, Quan = 0.95,
                     Plot = T, axisLimit = 40, dataPoint = 1)
jointOut$coverage

##  Find the bandwidth via cross-validation for the target rate of 0.95
target_rate=0.95
Thres_u=0.95
seqbw=seq(0.1,0.7,by=0.05) # a seq of bandwidth
cv_bw=rep(NA,length(seqbw))
## coverage rates for each target quantile
cv_result=sapply(seqbw,function(bw) crossValidate(Dat = rbind(fPrcp_tr,fPrcp_te),Ang =CP_prcp1$angular,
                                                  pMass = CP_prcp1$pmass,Thres = Thres_u,
                                                  kfold = 5,
                                                  bandW = bw,Quan = target_rate))
##  Find cv-bandwidths
cv_bw=seqbw[which.min(abs(cv_result-target_rate))]
cv_bw

####  Plot of KDE for angular densities
kde_out=KDE_w(Ang = CP_prcp1$angular,Pmass = CP_prcp1$pmass,bw = T,h = 0.3,Plot=T)
kde_out=KDE_w(Ang = CP_prcp1$angular,Pmass = CP_prcp1$pmass,bw = T,h = 0.4,Plot=T)
kde_out=KDE_w(Ang = CP_prcp1$angular,Pmass = CP_prcp1$pmass,bw = F,Plot=T) # bandwidth: 0.1181188
kde_out

####  Plot an approximate conditional density with the 95% conditional interval
dev.new()
conden_out=condDensity(xhatPoint=predVec[415],xp1Point=obsVec[415],
                       kde_h = kde_out, Quan = 0.95,
                       xlim_upper = 7,
                       ylim_upper = 10, Plot_h = FALSE)

####  Assess the coverage rate
####  Plot conditional intervals with lines to reproduce Figure 4 (right)
target_rate=0.95
XhatXp1 <- cbind(predVec,obsVec)
Keep <- XhatXp1[,1] > quantile(XhatXp1[,1], target_rate)
XhatXp1_top <- XhatXp1[Keep,]
par(mar=c(5.1,5.1,2,2))
coverOut=coverageRate(XY = XhatXp1_top, kde_est = kde_out, Quan = target_rate, Plot = F)
coverOut$CoverageRate # 0.84/0.95
dev.off()

