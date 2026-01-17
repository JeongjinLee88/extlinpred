##  Set your working directory
#setwd(...)
setwd("/home/leej40/Documents/extlinear/code")
source("functions.R")
library(Matrix) #nearPD
library(MASS) #ginv
library(plotrix)  # 'draw.circle'
library(VGAM) # probit
library(ks) # kde

##  Figure 1 in supplement: Calculate the tail ratio of prediction error
set.seed(1205)
PredictionError(Ncol = 10, Nrow = 4, n = 20000, min = 0, max = 5, plot_D=T)

##  Simulate from a shifted Pareto dist
set.seed(1234)
Nrow=7; Ncol=400; n=60000; min=0; max=5
U <- runif(n*Ncol)
shift <- 0.9352074  # make the mean of InvT(Z) centered.
Z <- matrix(1/sqrt(1-U)-shift,nrow=n,ncol=Ncol) # Necessary
##  Generate a p x q matrix A from a uniform dist
B <- matrix(runif(Nrow*Ncol,min = min, max = max), nrow = Nrow, ncol = Ncol)
B_norm=sqrt(apply(B^2,1,sum))
A <- B/B_norm

##  Generate X = A o Z
Out=genDataParams(A = A, Z = Z)
Data_X=Out$Xp
set.seed(12345)
R_index=sample(x = seq(1:n),size = n,replace = F)
Data_X=Data_X[R_index,]
##  Split the data set into a training / test set
n_train=ceiling(n*(2/3))
Train=Data_X[1:n_train,]
Test=Data_X[-(1:n_train),]

##  Estimate the TPDM, TPDM_pred, and b.
Thres_u=0.75  # tried 0.9, 0.95, 0.98
Est=estimateParams(X = Train, Thres = Thres_u)
Est$TPDM_hat
Out$TPDM_X # true TPDM

Est$TPDM_Phat
Out$TPDM_Pred # true prediction TPDM

##  CP-factorization for a 2x2 TPDM_Phat
CPout=CPfactor(Mtx = Est$TPDM_Phat,q_star = 10,ite_pereach = 5000,ite_cp = 100)
load(file="CPout.RData")

##  Create the 95% joint polar region from a CP-factorization
##  Calculate Xhat in the test set
Xhat_test=Amul(t(Est$bhat),t(Test[,-Nrow]))
Xhat_test=as.vector(Xhat_test)
##  Find the 95% joint polar region
target_rate=0.95
jointOut=jointRegion(Xhat = Xhat_test, Xf = Test[,Nrow],
                     Angular = CPout$angular, Pmass = CPout$pmass, Quan = target_rate,
                     Plot = T, axisLimit = 80, dataPoint = 14172)
jointOut$coverage
# 0.966
Xhat_test[14172] # 27.80396 (should be updated in the manuscript)
Test[14172,7] # 27.97971 (should be updated in the manuscript)

##  Kernel density estimation to approximate an angular density 'h'
##  Angular components and masses obtained from a CP-factor
kde_out=KDE_w(Ang = CPout$angular, Pmass = CPout$pmass, bw = F, Plot=T)
kde_out

##  Figure 2 (right): Plot an approximate conditional density with the 95% conditional interval
dev.new()
conden_out=condDensity(xhatPoint=Xhat_test[14172],xp1Point=Test[14172,Nrow],
                       kde_h = kde_out, Quan = target_rate,
                       xlim_upper = Xhat_test[14172]+0.06, ylim_upper = 90, Plot_h = T)

##  Assess the coverage rate
##  Figure 2 (left): Plot conditional intervals using the default bandwidth
target_rate=0.95
XhatXp1 <- cbind(Xhat_test,Test[,Nrow])
Keep <- XhatXp1[,1] > quantile(XhatXp1[,1], target_rate)
XhatXp1_top <- XhatXp1[Keep,]
coverOut=coverageRate(XY = XhatXp1_top, kde_est = kde_out, Quan = target_rate, Plot = T)
coverOut$CoverageRate
# 0.94


