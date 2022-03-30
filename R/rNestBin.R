#' Random nested exchangeable binary numbers 
#'
#' This function allows you to sample multivariate binary random variables with 
#'   a nested exchangeable correlation matrix, and different prevalences (means) 
#'   for each sub-cluster.
#'
#' The output data are binary random variables separated into clusters. Clusters 
#'   are independent, and each cluster has the same multivariate joint 
#'   distribution as all others. Each cluster is separated into C subclusters. 
#'   There are two levels of correlation: between subcluster and within 
#'   subcluster. All observations in the same subcluster have the same 
#'   prevalence, but different subclusters can have different prevalences.
#'   
#' The elements of the correlation matrix are given by 
#'   R[i,j] = 1 if i = j
#'          = rhoCT if i and j are observations in the same subcluster
#'          = rhoC if i and j are observations in the same cluster but not the 
#'                    same subcluster. rhoC <= rhoCT
#'                    
#' @param means Vector of prevalences, one for each subcluster within the cluster. The length of this vector determines the number of subclusters per cluster.
#' @param rhoC Correlation between two observations in the same cluster, but not the same subcluster.
#' @param rhoCT Correlation between two observations in the same subcluster.
#' @param n Number of observations per subcluster.
#' @param C Number of clusters to sample.
#' @param sample If sample=FALSE then just check whether the values of rhoC, rhoCT, means are feasible and stop if they are not feasible.
#' @export
#' @examples rNestBin(rhoC=0.2,rhoCT=0.4,means=rep(c(0.4,0.1),1),n=20,C=10,sample=T)
rNestBin <- function(means,rhoC,rhoCT,n=1,C=1,sample=TRUE){ 

  if(rhoc>rhoct){
    stop("rhoc must not be greater than rhoct")
  }
  if(rhoc<=0){
    stop("rhoc and rhoct must be positive")
  }
  
  TT = length(means) # number of subclusters (or time periods)
  NN = TT*n # total number of observations per cluster
  PPV = rep(means,each=n)
  
  # prentice constraints: maximum value for rhoC
  Pl = min(means)
  Pu = max(means)
  max_corr = sqrt(Pl*(1-Pu)/(Pu*(1-Pl)))
  if(rhoC>max_corr){
    stop('Prentice constraints violated')
  }
  
  # testing for feasibility
  lz = sqrt(rhoC*means*(1-means))
  ly = sqrt((rhoCT-rhoC)*means*(1-means))
  if(any(lz>0.5)){
    stop('between subcluster correlation rhoC is too high, or individual total variance is too high: mixing coefficient for Z is greater than 1')
  }
  if(any(ly>0.5)){
    stop('rhoCT-rhoC is too high, or individual total variance is too high: mixing coefficient for Y is greater than 1.')
  }
  if(any(lz+ly>0.5)){
    stop('sqrt(rhoCT-rhoC)+sqrt(rhoC) is too high, or individual total variance is too high: mixing coefficient for Y is greater than 1.')
  }
  
  # xi bounds section
  aa = sqrt(rhoCT-rhoC)
  bb = sqrt(rhoC)
  cc = sqrt(means/(1-means))
  det = cc^2*((aa^2-bb^2-1)^2-4*bb^2)
  if(any(det<0)){
    stop('No feasible values of z')
  }
  qzL = (-cc*(aa^2-bb^2-1)-sqrt(det))/(2*bb)
  qzU = (-cc*(aa^2-bb^2-1)+sqrt(det))/(2*bb)
  if(max(qzL)>min(qzU)){
    stop('No feasible values of z')
  }
  qz = (max(qzL)+min(qzU))/2 # sqrt odds for z
  qyL = aa/(1/cc-bb/qz)
  qyU = (cc-bb*qz)/aa
  qy = (qyL+qyU)/2# sqrt odds for y
  if(rhoc==rhoct){ # in this case aa=0, so the value of qy is irrelevant
    qy = rep(1,length(PP))
  }
  stopifnot(aa/qy+bb/qz<=1/cc)
  stopifnot(aa*qy+bb*qz<=cc)
  zz = qz^2/(1+qz^2)
  yy = qy^2/(1+qy^2)
  
  soy = sqrt(yy/(1-yy)) # sqrt odds for y, same as qy
  vary = yy*(1-yy) 
  soz = sqrt(zz/(1-zz)) # sqrt odds for z, same as qz
  varz = zz*(1-zz) 
  mz = lz/sqrt(varz) # length TT
  my = ly/sqrt(vary) # length TT
  # vector versions, length NN = n*TT
  yyV = rep(yy,each=n)
  zzV = rep(zz,each=NN)
  mzV = rep(mz,each=n)
  myV = rep(my,each=n)
  mxV = 1-mzV-myV
  soyV = rep(soy,each=n)
  varyV = rep(vary,each=n)
  sozV = rep(soz,each=NN)
  varzV = rep(varz,each=NN)
  xxV = (PPV-myV*yyV-mzV*zzV)/mxV
  
  if(sample){
    # sampling
    WW = matrix(0,nrow=C,ncol=NN)
    XX = matrix(0,nrow=C,ncol=NN)
    YY = matrix(0,nrow=C,ncol=NN)
    ZZ = matrix(0,nrow=C,ncol=NN)
    UU = matrix(0,nrow=C,ncol=NN)
    for(i in 1:C){
      XX[i,] = rbinom(n=NN,size=1,prob=xxV)
      YY[i,] = rep(rbinom(n=TT,size=1,prob=yy),each=n)
      ZZ[i,] = rep(rbinom(n=1,size=1,prob=zz),each=NN)
      UU[i,] = runif(n=NN,0,1)
      mixture = XX[i,]*(UU[i,]<=mxV)+YY[i,]*(UU[i,]>mxV)*(UU[i,]<=(mxV+myV))+ZZ[i,]*(UU[i,]>=(mxV+myV))
      WW[i,] = mixture
    }
  }else{
    WW=-1
  }
  return(WW)
}