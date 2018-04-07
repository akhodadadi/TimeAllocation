NLL_Model11 = function(x,parameter)
{
  #===parameters===
  #---free parameters---
  if (parameter$numOfThreshold==2){
    #     sigma01=x[1];sigma2=x[2];
    #     a01=x[3];a02=x[4];ap1=x[5];ap2=x[6];
    #     lambda1=x[7];lambda2=x[8];k1=x[9];k2=x[10];
    #     r1_p=x[11];r2_p=x[12];r1_n=x[13];r2_n=x[14];
    #     c1=x[15];c2=x[16];
  }else{
    sigma01=x[1];sigma02=sigma01;
    a01=x[2];a02=a01;ap1=x[3];ap2=ap1;
    lambda1=x[4];lambda2=lambda1;k1=x[5];k2=k1;
    r1_p=x[6];r2_p=r1_p;r1_n=x[7];r2_n=r1_n;
    c1=x[8];c2=c1;
    alpha_s1=x[9];alpha_s2=alpha_s1;
  }
  #---free parameters---
  
  #---other---
  rtVect=parameter$rtVect
  conditionVect=parameter$conditionVect;
  rewardVect=parameter$rewardVect;
  positionVect=parameter$positionVect;
  numOfTrials=length(rewardVect);
  whatToDo=parameter$whatToDo
  ind1=conditionVect==1;ind2=conditionVect==2;
  #---other---
  #===parameters===
  
  #===compute sigma===
  sigma1=sigma01-alpha_s1*seq(1,sum(ind1));
  sigma1[sigma1<5]=5;sigma1[sigma1>100]=100;
  sigma2=sigma02-alpha_s2*seq(1,sum(ind2));
  sigma2[sigma2<5]=5;sigma2[sigma2>100]=100;
  #===compute sigma===

  #===compute utility of reward===
  rewardVect[rewardVect>0 & ind1]=r1_p
  rewardVect[rewardVect>0 & ind2]=r2_p
  rewardVect[rewardVect<0 & ind1]=r1_n
  rewardVect[rewardVect<0 & ind2]=r2_n
  #===compute utility of reward===
    
  #===simulate model===
  aMat = matrix(0,nrow=2,ncol=numOfTrials);#'a' in weibull
  aMat[1,1]=a01;aMat[2,1]=a02;
  bMat=matrix(0,nrow=2,ncol=numOfTrials);#thrshold in each trial
  
  delta = matrix(0,nrow=1,ncol=numOfTrials);
  delta[ind1] = rewardVect[ind1] - rtVect[ind1]*c1;
  delta[ind2] = rewardVect[ind2] - rtVect[ind2]*c2;
  
  for (k in 1:(numOfTrials-1)){
    #---set parameteres for current condition---
    s=conditionVect[k]
    if (s==1){
      a=aMat[1,k];ap=ap1;lambda=lambda1;k_b=k1;#threshold params
      
    }else{
      a=aMat[2,k];ap=ap2;lambda=lambda2;k_b=k2;#threshold params
    }
    #---set parameteres for current condition---
    
    #---compute value of threshold---
    bMat[s,k] = a-(1-exp(-(rtVect[k]/lambda)^k_b))*(.5*a-ap);
    #---compute value of threshold---
    
    #---update 'a' in Weibull threhsolds---
    aMat[,k+1]=aMat[,k];
    aMat[s,k+1] = aMat[s,k] -
      delta[k]*(bMat[s,k]-positionVect[k])*(1+.5*(1-exp(-(rtVect[k]/lambda)^k_b)));
    
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
  bMat[s,k] = a-(1-exp(-(rtVect[k]/lambda)^k_b))*(.5*a-ap);
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
  }
  #===compute the likelihood===
  
}
