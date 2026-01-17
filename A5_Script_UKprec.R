##  Set a working directory
#setwd(...)
setwd("/home/leej40/Documents/extlinear/data")
load("cumbriaLanc.RData")
#should load four R objects: prcp, stations, dateVec, dat
#prcp is dat without first column of dates to make indexing more simple
prcp=as.data.frame(prcp)
stations=as.data.frame(stations)
dat=as.data.frame(dat)

##  Explore the data
numDays <- dim(prcp)[1]  # 56274
numStations <- dim(stations)[1]  # 608
pctMissing <- apply(is.na(prcp), 2, sum)/numDays
##  retaining data only after 1960
date1960 <- which(dateVec == "1960-01-01")  # 32506
prcp <- prcp[date1960:numDays,]
dat <- dat[date1960:numDays,]
numDays <- dim(prcp)[1]
pctMissing <- apply(is.na(prcp), 2, sum)/numDays

##  Retaining stations which have more missing < .1
keep <- pctMissing < .1
sum(keep)
prcp <- prcp[,keep]
dat <- dat[, c(1,which(keep)+1)]
stations <- stations[keep,]
numStations <- dim(stations)[1] # 30
pctMissing <- pctMissing[keep]
prcp=na.omit(prcp)
dim(prcp)
dat=na.omit(dat)
dim(dat)

##  Find max of each station
maxAtStat <- apply(prcp, 2, max, na.rm = T)
range(maxAtStat) # (49.4, 201.2)
##  look at day which is has absolute maximum as sanity check
statWithMax <- which(maxAtStat == max(maxAtStat)) # 26
dayOfMax <- which(prcp[, statWithMax] == max(maxAtStat)) #3375
maxDayPrcp <- prcp[dayOfMax,]
library(fields)
# Convert the 1-row tibble to a numeric vector
vals <- as.numeric(maxDayPrcp)
# Convert the coordinate tibble to a matrix
coords <- as.matrix(stations[, c("lon", "lat")])
# Run as.image with the correct data types
ai <- as.image(vals, x = coords)
dev.new()
image.plot(ai)

##  Shape xi estimates from GPD
xiEst <- numeric(numStations)
library(lmomco)
for(i in seq(1, numStations)){
  prcpi <- prcp[!is.na(prcp[,i]),i]
  q.99 <- quantile(prcpi, .99)
  prcpiExc <- prcpi[prcpi > q.99]
  lmomEsts <- lmoms(prcpiExc)
  xiEst[i] <- -pargpa(lmomEsts)$para[3]
}  
mean(xiEst) # 0.07
mean(xiEst>0)  # 80% have positive estimates
dev.new()
plot(xiEst)  
abline(h = 0, lty = 3)

##  Divide the data set into two sets randomly
##  Rank transform to Pareto-2 marginals
set.seed(12345)
n_obs=dim(prcp)[1]
R_index=sample(x = seq(1:n_obs),size = n_obs,replace = F)
prcp=prcp[R_index,]
trun=ceiling(n_obs*(2/3))
prcp_tr=prcp[1:trun,] # training set
prcp_te=prcp[-(1:trun),] # test set

U_tr <- apply(prcp_tr, 2, function(i) rank(i)/(dim(prcp_tr)[1]+1))
U_te <- apply(prcp_te, 2, function(i) rank(i)/(dim(prcp_te)[1]+1))

##  Transformed observations approximating a shifted Pareto dist
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
    rad_thres <- quantile(rad, 0.95)
    radKeep <- rad > rad_thres
    angKeep <- ang[radKeep,]
    tpdEst <- 2*mean(apply(angKeep, 1, prod))
    TPDM[i,j] <- TPDM[j,i] <- tpdEst
  }
}
diag(TPDM) <- 1

##  Figure 4 (left): Map showing the strength of dependence (by TPDM)
##  with location 2, chosen as it has lowest pctMissing
loi <- which(pctMissing == min(pctMissing))  # coniston-garden-house (labelled 2)
ai <- as.image(TPDM[loi,], x = coords)
pdf("/home/leej40/Documents/extlinear/TPDM_final.pdf",6,6)
dev.new()
par(mar = c(5, 5, 4, 1))
image.plot(ai,cex.main=1.5, cex.lab=1.5, cex.axis=1.5, cex=1.5,axis.args = list(cex.axis = 1.2),xlab="lon",ylab="lat")
points(stations[loi, "lon"]+0.01, stations[loi, "lat"],pch=8, cex=1,lwd=1)
dev.off()

##  Figure 4 (right): Explore optimal weights b*
idx <- 1:numStations
idx <- c(idx[-loi], idx[loi])
TPDM_reord <- TPDM[idx, idx]
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

##  Load necessary functions for prediction
setwd("/home/leej40/Documents/extlinear/code")
source("functions.R")
library(MASS) #ginv
library(plotrix)  # 'draw.circle'
library(VGAM)
library(ks)

##  Perform prediction in a test set
numPred <- dim(fPrcp_te)[1]
predVec <- numeric(numPred)
obsVec <- numeric(numPred)

for(i in seq(1:numPred)){
  idx_loi <- as.numeric(loi)
  P_reduced <- dim(TPDM)[1]
  idx <- seq(1, P_reduced)
  idx <- c(idx[-idx_loi], idx_loi)
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

##  Reorder the indices of the TPDM
idx <- seq(1, dim(TPDM)[1])
idx <- c(idx[-loi], loi)
TPDM_reord <- TPDM[idx, idx]
d=dim(TPDM_reord)[1]
Sig11=TPDM_reord[1:(d-1),1:(d-1)]
Sig12=TPDM_reord[1:(d-1),d]
Sig21=t(Sig12)
Sig22=TPDM_reord[d,d]
##  Prediction inner product matrix
Sig_P_11=Sig21%*%solve(Sig11)%*%Sig12
TPDM_Phat=matrix(c(Sig_P_11,Sig_P_11,Sig_P_11,Sig22),nrow=2)

##  Reiterate the decomposition of the TPDM_Pred_hat (it takes a bit time)
##  To skip it, load the saved RData.
CP_prcp1=CPfactor(Mtx = TPDM_Phat,q_star = 10,ite_pereach = 5000,ite_cp = 100)
save(CP_prcp1,file="CP_prcp_final.RData")
load(file="CP_prcp_final.RData")

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

##  Plot of KDE for angular densities
kde_out=KDE_w(Ang = CP_prcp1$angular,Pmass = CP_prcp1$pmass,bw = T,h = 0.4,Plot=T)
kde_out=KDE_w(Ang = CP_prcp1$angular,Pmass = CP_prcp1$pmass,bw = F,Plot=T) # bandwidth: 0.12
kde_out

##  Plot an approximate conditional density with the 95% conditional interval
dev.new()
conden_out=condDensity(xhatPoint=predVec[20],xp1Point=obsVec[20],
                       kde_h = kde_out, Quan = 0.95,
                       xlim_upper = 7,
                       ylim_upper = 10, Plot_h = FALSE)

##  Assess the coverage rate
target_rate=0.95
XhatXp1 <- cbind(predVec,obsVec)
Keep <- XhatXp1[,1] > quantile(XhatXp1[,1], target_rate)
XhatXp1_top <- XhatXp1[Keep,]
par(mar=c(5.1,5.1,2,2))
coverOut=coverageRate(XY = XhatXp1_top, kde_est = kde_out, Quan = target_rate, Plot = F)
coverOut$CoverageRate # 0.84/0.955
dev.off()

