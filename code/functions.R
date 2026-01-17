##  1. Tranformed linear operations
##  Mapping reals to positive reals
Trans = function(y){
  x=log(exp(y)+1)
  x[!is.finite(x)]=y[!is.finite(x)]
  return(x)
}
##  Mapping positive reals to reals
InvT = function(x){
  y=log(exp(x)-1)
  y[!is.finite(y) & x > 1]=x[!is.finite(y) & x > 1]
  return(y)
}
##  Transformed matrix multiplication
Amul <- function(A, x){
  y=InvT(x)
  Ay=A%*%y
  Ay=Trans(Ay)
  return(Ay)
}
##  Transformed scalar multiplication
Cmul <- function(c, x){
  y <- InvT(x)
  cy <- c*y
  cx <- Trans(cy)
  return(cx)
}
##  Transformed addition
vSum <- function(v1, v2){
  y1 <- InvT(v1)
  y2 <- InvT(v2)
  sumX <- y1 + y2
  sumV <- Trans(sumX)
  return(sumV)
}
##  Transformed subtraction
vSub <- function(v1, v2){
  y1 <- InvT(v1)
  y2 <- InvT(v2)
  sumX <- y1 - y2
  sumV <- Trans(sumX)
  return(sumV)
}

##  2. Generate data via a matrix multiplication 
PredictionError <- function(Ncol=10,Nrow=4,n=20000,min=0,max=5,plot_D=T){
  
  ##  Simulate a shifted Pareto dist
  U <- runif(n*Ncol)
  shift <- 0.9352074  # make the mean of InvT(Z) centered.
  Z <- matrix(1/sqrt(1-U)-shift,nrow=n,ncol=Ncol)
  
  ##  Generate a standardized matrix A from a uniform dist
  B <- matrix(runif(Nrow*Ncol,min = min, max = max), nrow = Nrow, ncol = Ncol)
  B_norm=sqrt(apply(B^2,1,sum))
  A <- B/B_norm
  
  ##  A known TPDM
  TPDM_X=A%*%t(A)
  
  ##  Tail ratio of prediction errors
  kstar=TPDM_X[Nrow,Nrow]-t(TPDM_X[1,2:Nrow])%*%solve(TPDM_X[1:(Nrow-1),1:(Nrow-1)])%*%TPDM_X[1,2:Nrow]
  
  ##  Simulate a random vector X by matrix multiplication
  X_t <- t(Amul(A, t(Z)))
  
  ##  Find the optimized vector b
  b=(solve(A[1:(Nrow-1),]%*%t(A[1:(Nrow-1),]),tol = 1e-30)%*%A[1:(Nrow-1),]%*%A[Nrow,]) # A[Nrow,] is a qx1 vector, not a 1xq.
  
  ##  The best transformed-linear predictor using TPDM
  Xhat=Amul(t(b),t(X_t[,1:(Nrow-1)]))
  Xhat=as.vector(Xhat)
  
  ##  Find the transformed-prediction error between Xhat and X_4
  Diff1=vSub(v1 = X_t[,Nrow],v2 = Xhat)
  Diff2=vSub(v1 = Xhat,v2 = X_t[,Nrow])
  D_TPDM=apply(X = cbind(Diff1,Diff2),1,max)
  
  ##  A scatterplot of "D" against Xhat with the 0.95 quantile for "D"
  if(plot_D){
    dev.new()
    par(mar=c(5.1,5.1,2,2))
    plot(Xhat,D_TPDM,xlim=c(0,max(Xhat)),ylim=c(0,max(Xhat)),main="",xlab=expression(hat(X)[4]),ylab="D", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    #plot(Xhat,D_TPDM,xlim=c(0,max(Xhat)),ylim=c(0,65),main="",xlab=expression(hat(X)[4]),ylab="D", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(h=sqrt(2*kstar/0.05),lty=2,lwd=2)
    abline(h=sqrt(2*kstar/0.05),lty=2,lwd=2)
  }
  
  ##  Assess the coverage rate
  D_up=as.numeric(sqrt(kstar/0.025))
  Coverage_D=1-2*length(D_TPDM[D_TPDM>D_up])/length(D_TPDM)
  
  return(list("kstar"=kstar,"b"=b,"D_TPDM"=D_TPDM,"Coverage_D"=Coverage_D))
}

##  3. Prediction (Xhat) / vector b / prediction TPDM / angular mass
genDataParams <- function(A, Z){
  
  p=dim(A)[1]
  ##  Specified TPDM
  TPDM_X=A%*%t(A)
  
  ##  Generate a vector X by matrix multiplication
  Xp <- t(Amul(A, t(Z)))
  
  ##  Optimal vector b = Sig11^(-1) * Sig12
  b=solve(TPDM_X[1:(p-1),1:(p-1)])%*%TPDM_X[1:(p-1),p] 
  
  ##  The best transformed-linear predictor
  Xhat=Amul(t(b),t(Xp[,1:(p-1)]))
  Xhat=as.vector(Xhat)
  
  ##  The predictor TPDM of Xhat and X
  A_Pred=rbind(t(b)%*%A[1:(p-1),],A[p,])  #2xq
  TPDM_Pred=A_Pred%*%t(A_Pred)
  
  ##  The squared scale of D
  PredE=TPDM_X[p,p]-TPDM_X[p,1:(p-1)]%*%solve(TPDM_X[1:(p-1),1:(p-1)])%*%TPDM_X[p,1:(p-1)]
  PredE
  
  ##  Known angular components and masses
  Ang_T=atan(A_Pred[2,]/A_Pred[1,])
  Ang_mass=apply(A_Pred^2,2,sum)
  
  return(list("Xp"=Xp, "TPDM_X"=TPDM_X, "Xhat"=Xhat, "b"=b, "A_Pred"=A_Pred, "TPDM_Pred"=TPDM_Pred, "PredE"=PredE, "Ang_T"=Ang_T, "Ang_mass"=Ang_mass))  
}

##  4. Pairwise TPDM estimation
estimateParams <- function(X,Thres=0.95){
  
  d=dim(X)[2]  
  
  ## Pairwise estimation
  TPDM_hat=matrix(0,d,d)
  for(i in 1:d){
    for(j in 1:i){
      R_pair=sqrt(apply(X[,c(i,j)]^2,1,sum))
      W_pair=X[,c(i,j)]/R_pair
      Keep=R_pair > quantile(R_pair, Thres)
      W_keep=W_pair[Keep,]
      TPDM_hat[i,j]=t(W_keep[,1])%*%W_keep[,2] / dim(W_keep)[1] * 2
    }
  }
  ##  Estimated TPDM
  TPDM_hat=TPDM_hat+t(TPDM_hat)-diag(diag(TPDM_hat))
  
  Sig11=TPDM_hat[1:(d-1),1:(d-1)]
  Sig12=TPDM_hat[1:(d-1),d]
  Sig21=t(Sig12)
  Sig22=TPDM_hat[d,d]
  
  ##  Prediction inner product matrix
  Sig_P_11=Sig21%*%solve(Sig11)%*%Sig12
  TPDM_Phat=matrix(c(Sig_P_11,Sig_P_11,Sig_P_11,Sig22),nrow=2)
  ##  b hat
  bhat=solve(Sig11)%*%Sig12
  
  return(list("TPDM_hat"=TPDM_hat,"TPDM_Phat"=TPDM_Phat,"bhat"=bhat))
}

##  5. Completely positive decomposition
CPfactor <- function(Mtx, q_star=9, ite_pereach=5000, ite_cp=80){
  
  dcol=dim(Mtx)[2]
  if(q_star < dcol){
    stop("q_star must be greater than  or equal to the column dimension of the given matrix")
  }
  add_col=q_star-(dcol-1)
  
  ##  TPDM_Pred (Rank=2)
  if(any(eigen(x = Mtx)$values < 0)){
    Mtx=nearPD(Mtx)
  }
  
  ##  Cholesky factorization for TPDM_Pred=BB' (Initial factorization)
  Btilde=t(chol(Mtx)) 
  #t(chol(x = Mtx))%*%chol(x = Mtx) = TPDM_Pred
  
  ##  Add extra columns to B
  AugmentedMatrix=matrix(rep(Btilde[,1],add_col)/sqrt(add_col),nrow=dim(Btilde)[1])
  B=cbind(Btilde[,2],AugmentedMatrix)
  #B%*%t(B) = TPDM_Pred
  
  ##  CP-factorization
  BQ_save=array(0,dim=c(dcol,q_star,ite_cp))
  for(j in 1:ite_cp){
    converged <- FALSE
    while(!converged){
      ##  Generate a qxq random matrix M
      M=matrix(rnorm((q_star)*(q_star)),nrow=(q_star))
      ##  Find Q_0 by projecting M onto O_r
      U=eigen(M%*%t(M))$vectors
      V=eigen(t(M)%*%M)$vectors
      ##  Initial orthogonal matrix
      Q_0=U%*%t(V)
      
      ##  Alternating projections
      BQ=B%*%Q_0
      for(i in 1:ite_pereach){
        BQ[BQ<=0]<-0
        D=BQ
        PseInv=ginv(X = D)
        Phat=PseInv%*%D+(diag((q_star))-PseInv%*%B)%*%Q_0
        LeftSingular=eigen(Phat%*%t(Phat))$vectors
        RightSingular=eigen(t(Phat)%*%Phat)$vectors
        Q_0=LeftSingular%*%t(RightSingular)
        BQ=B%*%Q_0
        BQ=Re(BQ)
        if(all(BQ >= -10^(-15))){
          converged <- TRUE
          break
        }
        print(i)
      }
    }
    BQ_save[,,j]=BQ
  }
  
  array_reorder=aperm(BQ_save,c(1,3,2)) # reorder dimensions to (dcol x ite_cp x q_star)
  BQ_comb=matrix(array_reorder, nrow=dcol)
  
  ##  save angular points and point masses
  pmass=apply(BQ_comb^2,2,sum)/ite_cp
  angular=BQ_comb[1,]/sqrt(apply(BQ_comb^2,2,sum))
  
  return(list("angular"=angular,"pmass"=pmass,"BQ_comb"=BQ_comb))  
}

##  6. Joint region
jointRegion <- function(Xhat, Xf, Angular, Pmass, Quan=0.95, Plot=TRUE, axisLimit=40, dataPoint){
  
  ##  Angular measure from CP-factorization
  Ang_c=cbind(Angular,Pmass)
  sl=sort.list(Angular)
  Ang_sort=Ang_c[sl,]
  Cum_s=cumsum(Ang_sort[,2])
  Ang_H=cbind(Ang_sort,Cum_s)
  
  ##  Total mass
  Tmass=sum(Pmass)
  
  ##  Find lower and upper quantiles of angular
  Lquan=(1-Quan)/2
  Uquan=1-Lquan
  wQuan=approx(y=Ang_H[,1], x=Ang_H[,3], xout=c(Tmass*Lquan,Tmass*Uquan), method = "linear", rule=2)$y
  
  ##  Transforming (Xhat, Xf) to polar coordinates
  XY <- cbind(Xhat,Xf)
  Rad <- sqrt(apply(XY^2, 1, sum))
  W <- XY/Rad
  Keep <- Rad > quantile(Rad, Quan)
  W_top <- W[Keep,]
  XY_top <- XY[Keep,]
  
  interiorRegion=W_top[W_top[,1]>=wQuan[1] & W_top[,1]<=wQuan[2], 1]
  coverage=length(interiorRegion)/dim(W_top)[1]
  
  ##  Plot of (Xhat, Xf) with a joint polar region
  X_comb=data.frame(XY,Rad)
  X_comb$color=as.character(cut(X_comb[,3],breaks = c(0,as.numeric(quantile(X_comb[,3],Quan)),Inf),labels = c("grey","black"),right = FALSE))
  R_thres=quantile(Rad,probs = Quan)  
  
  if(Plot){
    dev.new()
    #pdf("/home/leej40/Documents/extlinear/simJointRegion.pdf",6,6)
    par(mar=c(5.1,5.1,2,2))
    plot(X_comb[,1],X_comb[,2],col=X_comb$color,main="",xaxs='i',yaxs='i',
         xlim=c(0,axisLimit),ylim=c(0,axisLimit),
         xlab=expression(paste("Predicted"," ",NO[2]," ","(Pareto Scale)")),
         ylab=expression(paste("Observed"," ",NO[2]," ","(Pareto Scale)")),
         cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1)
    X_comb=X_comb[-dataPoint,]
    points(x = XY[dataPoint,1],y=XY[dataPoint,2],pch="*",col="blue",cex=4,lwd=2)
    draw.circle(x = 0, y = 0, radius = R_thres, border = "red", lwd=1, lty=1)
    abline(a = 0, b = sqrt(1-wQuan[1]^2)/wQuan[1], col="blue", lwd=2, lty=2)
    abline(a = 0, b = sqrt(1-wQuan[2]^2)/wQuan[2], col="blue", lwd=2, lty=2)
    #dev.off()
  }
  
  return(list("Ang_H"=Ang_H,"XY"=XY,"XY_top"=XY_top,"Tmass"=Tmass,"wQuan"=wQuan,"coverage"=coverage))
}

##  7. kernel density estimate
KDE_w <- function(Ang, Pmass, bw=FALSE, h=0.1, Plot=TRUE){
  
  wseq=seq(0,1,len=5000)
  
  ##  Transformed kde for h(w) at wseq
  ##  probitlink: use 'bvalue' to prevent numerical instability at close to 0 or 1.
  suppressWarnings(
    if(bw){
      kde_trans=ks::kde(x = probitlink(Ang,bvalue = .Machine$double.eps), w = Pmass, eval.points = probitlink(wseq,bvalue = .Machine$double.eps), h=h)
    }else{
      kde_trans=ks::kde(x = probitlink(Ang,bvalue = .Machine$double.eps), w = Pmass, eval.points = probitlink(wseq,bvalue = .Machine$double.eps))  
    }
  )
  
  ##  Adjust boundary bias
  kde_trans$eval.points=probitlink(kde_trans$eval.points, bvalue = .Machine$double.eps, inverse = T)
  kde_trans$estimate=(kde_trans$estimate)*(1/dnorm(kde_trans$eval.points))
  
  ##  Plot a kernel density in terms of w_L2
  if(Plot){
    dev.new()
    par(mar=c(5.1,5.1,2,2))
    plot(kde_trans,xlim=c(0,1),main="",xlab="w",ylab="",
         cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
    points(Ang,Pmass)
  }
  
  kde_w=kde_trans$estimate
  return(kde_w)
}

##  8. Cross-validated bandwidth for kernel density
crossValidate <- function(Dat, Ang, pMass, Thres, kfold=3, bandW, Quan){
  
  n=nrow(Dat);p=dim(Dat)[2]
  fold_size=floor(n/kfold)
  
  coverRate=c()
  
  for(i in 1:kfold){
    test_ind=((i-1)*fold_size+1):(i*fold_size)
    if(i==kfold){
      test_ind=((i-1)*fold_size+1):n
    }
    
    Train=Dat[-test_ind,]
    Test=Dat[test_ind,]
    
    Est=estimateParams(X = Train, Thres = Thres)
    b_vec=Est$bhat
    
    Xhat_test=Amul(t(b_vec),t(Test[,-p]))
    Xhat_test=as.vector(Xhat_test)
    jointOut=jointRegion(Xhat = Xhat_test, Xf = Test[,p],
                         Angular = Ang, Pmass = pMass, Quan = Quan,
                         Plot = F, axisLimit = 40, dataPoint = 1)
    
    if(bandW==0){
      kde_out=KDE_w(Ang = Ang,Pmass = pMass,Plot=F)
    }else{
      kde_out=KDE_w(Ang = Ang,Pmass = pMass,bw = T,h = bandW,Plot=F)
    }
    #coverOut=coverageRate(XY = jointOut$XY_top, kde_est = kde_out, Quan = Quan,Plot = F)
    
    XhatXp1 <- cbind(Xhat_test,Test[,p])
    Keep <- XhatXp1[,1] > quantile(XhatXp1[,1], Quan)
    XhatXp1_top <- XhatXp1[Keep,]
    coverOut=coverageRate(XY = XhatXp1_top, kde_est = kde_out, Quan = Quan,Plot = F)
    coverRate<-c(coverRate,coverOut$CoverageRate)
  }
  ##  Average coverage rate across folds
  mean(coverRate)
}

##  9. Prediction interval given a large predictand
condDensity <- function(xhatPoint, xp1Point, kde_h, Quan, xlim_upper, ylim_upper, Plot_v = TRUE, Plot_h = TRUE){
  
  ##  Calculate the approximate conditional density for a given large Xhat
  ##  Define angular components 
  wseq=seq(1,1e-10,len=length(kde_h)) 
  
  ##  Radial components and corresponding x_p+1
  rad <- xhatPoint/wseq 
  xp1 <- rad*(1-wseq) 
  
  #xp1<-xhatPoint*sqrt(wseq^2/(1-wseq^2))
  #rad<- sqrt(xhatPoint^2+xp1^2)
  
  ##  'denApproxNorm' returns an approximated angular density from CPD
  denApproxNoNorm=2*(rad^(-5))*rev(kde_h)*xp1 
  denApproxNoNorm[is.nan(denApproxNoNorm)] <- 0
  
  xp1Temp <- xp1
  #xp1Temp[length(xp1Temp)] <- 2 * xp1Temp[length(xp1Temp) - 1]
  delta <- xp1Temp[2:(length(xp1Temp))] - xp1Temp[1:(length(xp1Temp) - 1)]
  piece1 <- delta * denApproxNoNorm[1:(length(denApproxNoNorm) - 1)]
  piece2 <- delta * denApproxNoNorm[2:(length(denApproxNoNorm))]
  traps <- 1/2*(piece1 + piece2)
  normalizer <- sum(traps)
  cumTraps <- c(0, cumsum(traps))/normalizer # CDF
  denApproxNorm <- denApproxNoNorm/normalizer
  
  
  ##  Sort the CDF
  Con_comb=cbind(xp1,cumTraps)
  sl=sort.list(xp1)
  Con_sort=Con_comb[sl,]
  ##  Interpolate the empirical CDF to find the 95% conditional interval
  Lquan=(1-Quan)/2;Uquan=1-Lquan
  condQuan=approx(x = Con_sort[,2], y = Con_sort[,1], xout = c(Lquan,Uquan), method = "linear",rule = 2,ties = "mean")$y
  #condQuan=approxExtrap(x = Con_sort[,2], y = Con_sort[,1], xout = c(Lquan,Uquan), method = "linear",rule = 2,ties = "mean")$y
  
  ##  Counting
  count=ifelse(xp1Point>=condQuan[1] & xp1Point<=condQuan[2],1,0)
  
  ##  Plot of the approximate conditional density with the 95% CI
  
  if(Plot_h){
    #dev.new()
    par(mar=c(5.1,5.1,2,2))
    plot(xp1, denApproxNorm, type = 'l', 
         ylim = range(denApproxNorm), xlim=c(0,50),
         xlab=expression(X[p+1]~"|"~hat(X)[p+1]),
         ylab = "Conditional density",
         main=expression(paste("The"," ","approximated"," ","conditional"," ","density"," ","when"," ",hat(X)[p+1]," ","is"," ","large")))
  }
  if(Plot_v){
    #dev.new()
    #pdf("/home/leej40/Documents/extlinear/simCondIntvlXHat27.pdf",6,6)
    par(mar=c(5.1,5.1,2,2))
    plot(xhatPoint+denApproxNorm, xp1, type='l', lwd=2,
         xlim=c(xhatPoint-0.01,xlim_upper),ylim=c(0,ylim_upper),
         xlab="Predicted values", ylab="Observed values",
         cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1)
    abline(v = xhatPoint,lwd=2)
    points(x = xhatPoint, y=xp1Point, pch="*", col="blue", cex=3, lwd=2)
    ##  95% conditional interval  
    points(x = xhatPoint,y=condQuan[1],pch="-", col="red",cex=3,lwd=2)
    points(x = xhatPoint,y=condQuan[2],pch="-", col="red",cex=3,lwd=2)
    #dev.off()
  }
  
  return(list("denApproxNoNorm"=denApproxNoNorm, "denApproxNorm"=denApproxNorm, "xp1"=xp1, "cumTraps"=cumTraps, "condQuan"=condQuan, "count"=count))
}

##  10. Calculate coverage rate
coverageRate <- function(XY, kde_est, Quan, Plot=TRUE){
  
  XY_top=XY
  
  ##  Assess the coverage rate
  counting=rep(0,dim(XY_top)[1])
  PI=matrix(0,nrow = dim(XY_top)[1],ncol = 2,byrow = T)
  for(i in 1:dim(XY_top)[1]){
    out=condDensity(xhatPoint=XY_top[i,1],xp1Point=XY_top[i,2],
                    kde_h = kde_est, Quan = Quan,
                    Plot_v = FALSE, Plot_h = FALSE)
    counting[i]=out$count
    PI[i,]=out$condQuan
    print(i)
  }
  CoverageRate=sum(counting)/dim(XY_top)[1]
  
  ##  Plot of (Xhat,X_p+1)
  ABC=cbind(XY_top[,1],PI)
  sl <- sort.list(ABC[,1])
  ABCSort=ABC[sl,]
  if(Plot){
    dev.new()
    pdf("/home/leej40/Documents/extlinear/PrcpParetoScale.pdf",6,6)
    #pdf("/home/leej40/Documents/extlinear/simCondIntvls.pdf",6,6)
    par(mar=c(5.1,5.1,2,2))
    plot(XY_top[,1], XY_top[,2], xlim=c(0,50), ylim=c(0,50),
         main="", xlab="", ylab=expression(X[p+1]),
         cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1.3)
    title(xlab=expression(hat(X)[p+1]),line=3.5, cex.lab=1.5)
    lines(ABCSort[,1],ABCSort[,2], lty = 2, cex=2, lwd=2, col="blue")
    lines(ABCSort[,1],ABCSort[,3], lty = 2, cex=2, lwd=2, col="blue")
    #points(x = XY_top[,1],y=PI[,1],pch="-", col="blue",cex=2,lwd=2)
    #points(x = XY_top[,1],y=PI[,2],pch="-", col="blue",cex=2,lwd=2)
    dev.off()
  }
  
  return(list("CoverageRate"=CoverageRate,"PI"=PI))
}




