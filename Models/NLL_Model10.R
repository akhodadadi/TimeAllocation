NLL_Model10 = function(x,parameter)
{
  #===parameters===
  #---free parameters---
  sigma01=x[1];sigma02=x[2];
  a01=x[3];a02=x[4];ap1=x[5];ap2=x[6];
  lambda1=x[7];lambda2=x[8];k1=x[9];k2=x[10];
  
  alpha_c=x[11];alpha_m1=x[12];alpha_m2=x[13];
  c1=x[14];c2=x[15];
  V10=x[16];V20=x[17];
  
  alpha_s1=x[18];alpha_s2=x[19];
  rhoWin = floor(x[20]);#window size for estimating rho
  #---free parameters---
  
  #---other---
  conditionVect=parameter$conditionVect;
  ind1=conditionVect==1;ind2=conditionVect==2;
  ind=ind1 | ind2
  conditionVect=conditionVect[ind];
  rewardVect=parameter$rewardVect[ind];
  positionVect=parameter$positionVect[ind];
  rtVect=parameter$rtVect[ind];
  ind1=ind1[ind];ind2=ind2[ind];
  numOfTrials=length(rewardVect);
  whatToDo=parameter$whatToDo
  #---other---
  #===parameters===
  
  #===compute sigma===
  sigma1=sigma01-alpha_s1*seq(1,sum(ind1));
  sigma1[sigma1<5]=5;sigma1[sigma1>100]=100;
  sigma2=sigma02-alpha_s2*seq(1,sum(ind2));
  sigma2[sigma2<5]=5;sigma2[sigma2>100]=100;
  #===compute sigma===
  
  #===simulate model===
  bMat = matrix(0,nrow=2,ncol=numOfTrials);#thrshold in each trial
  aMat = matrix(0,nrow=2,ncol=numOfTrials);#'a' in weibull
  aMat[1,1]=a01;aMat[2,1]=a02;
  rho=rep(0,numOfTrials-1);#estimated average reward
  
  vMat = matrix(0,nrow=2,ncol=numOfTrials);
  vMat[1,1]=V10;vMat[2,1]=V20;
  
  for (k in 1:(numOfTrials-1)){
    
    #---estimate reward rate---
    if (k>rhoWin){
      rho[k]=
        sum(rewardVect[(k-rhoWin):k]) / sum(rtVect[(k-rhoWin):k])
    }else if (k<=rhoWin){
      rho[k]=
        sum(rewardVect[1:k]) / sum(rtVect[1:k])
    }
    #---estimate reward rate---
    
    #---states at k and k+1---
    s=conditionVect[k];sPrime=conditionVect[k+1];
    #---states at k and k+1---
    
    #---set parameteres for current condition---
    if (s==1){
      alpha_m=alpha_m1;c=c1;
      a=aMat[1,k];ap=ap1;lambda=lambda1;k_b=k1;#threshold params
    }
    else{
      alpha_m=alpha_m2;c=c2;
      a=aMat[2,k];ap=ap2;lambda=lambda2;k_b=k2;#threshold params
    }
    #---set parameteres for current condition---
    
    #---compute delta---
    delta = rewardVect[k] - c*rtVect[k]*rho[k] +
      vMat[sPrime,k] - vMat[s,k];
    #---compute delta---
    
    #---update state values---
    vMat[,k+1]=vMat[,k];
    vMat[s,k+1] = vMat[s,k] + alpha_c * delta;
    #---update state values---
    
    #---compute value of threshold---
    bMat[s,k] = a-(1-exp(-(rtVect[k]/lambda)^k_b))*(.5*a-ap);
    #---compute value of threshold---
    
    #---update 'a' in Weibull threhsolds---
    aMat[,k+1]=aMat[,k];
    aMat[s,k+1] = aMat[s,k] - alpha_m * delta * 
      (bMat[s,k]-positionVect[k])*(1+exp(-(rtVect[k]/lambda)^k_b));
    
    if (aMat[s,k+1]<0){
      aMat[s,k+1]=0;
    }
    if (aMat[s,k+1]>300){
      aMat[s,k+1]=300;
    }
    #---update 'a' in Weibull threhsolds---
  }
  
  #---compute threshold for last trial---
  k=numOfTrials
  s=conditionVect[k]
  if (s==1){
    a=aMat[1,k];ap=ap1;lambda=lambda1;k_b=k1;#threshold params
  }
  else{
    a=aMat[2,k];ap=ap2;lambda=lambda2;k_b=k2;#threshold params
  }
  bMat[s,k] = 
    a-(1-exp(-(rtVect[k]/lambda)^k_b))*(.5*a-ap);
  #---compute threshold for last trial---
  #===simulate model===
  
  #===compute the likelihood===
  if (whatToDo==1){#fitting
    n1=sum(ind1);n2=sum(ind2);
    SSE1= sum( ((positionVect[ind1]-bMat[1,ind1])^2)/(2*sigma1^2) );
    SSE2= sum( ((positionVect[ind2]-bMat[2,ind2])^2)/(2*sigma2^2) );
    
    NLL1=sum(log(sigma1))  + SSE1;
    NLL2=sum(log(sigma2))  + SSE2;
    NLL = NLL1+NLL2+.5*numOfTrials*log(2*pi);
  }else if (whatToDo==2){#simulate fitted
    output = list(b1Vect=bMat[1,],b2Vect=bMat[2,])
  }else if (whatToDo==3){#output param. of threshold
    threshold_param_mat=cbind(aMat[1,],
                              rep(ap1,numOfTrials),
                              rep(lambda1,numOfTrials),
                              rep(k1,numOfTrials),
                              aMat[2,],
                              rep(ap2,numOfTrials),
                              rep(lambda2,numOfTrials),
                              rep(k2,numOfTrials)
    )#this matrix keeps 
    #parameters of Weibull threshold in each trial.
    rownames(threshold_param_mat)=1:nrow(threshold_param_mat)
    colnames(threshold_param_mat)=c('a1','ap1','lambda1','k1',
                                    'a2','ap2','lambda2','k2')
    output=list(threshold_param_mat=threshold_param_mat)
  }
  #===compute the likelihood===
  
}
