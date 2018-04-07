NLL_Model8 = function(x,parameter)
{
  #===parameters===
  #---free parameters---
  sigma01=x[1];sigma02=x[2];
  m10=x[3];m20=x[4];
  alpha_c=x[5];alpha_m1=x[6];alpha_m2=x[7];
  c1=x[8];c2=x[9];
  V10=x[10];V20=x[11];
  
  alpha_s1=x[12];alpha_s2=x[13];
  #---free parameters---
  
  #---other---
  conditionVect=parameter$conditionVect;
  ind1=conditionVect==1;ind2=conditionVect==2;
  ind=ind1 | ind2
  conditionVect=conditionVect[ind];
  rewardVect=parameter$rewardVect[ind];
  positionVect=parameter$positionVect[ind];
  rtVect=parameter$rtVect[ind]
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
  mMat = matrix(0,nrow=2,ncol=numOfTrials);
  mMat[1,1]=m10;mMat[2,1]=m20;
  
  vMat = matrix(0,nrow=2,ncol=numOfTrials);
  vMat[1,1]=V10;vMat[2,1]=V20;
  
  for (k in 1:(numOfTrials-1)){
    #---states at k and k+1---
    s=conditionVect[k];sPrime=conditionVect[k+1];
    #---states at k and k+1---
    
    #---set parameteres for current condition---
    if (s==1){
      alpha_m=alpha_m1;c=c1;
    }
    else{
      alpha_m=alpha_m2;c=c2;
    }
    #---set parameteres for current condition---
    
    #---compute delta---
    delta = rewardVect[k] - c*rtVect[k] +
      vMat[sPrime,k] - vMat[s,k];
    #---compute delta---
    
    #---update state values---
    vMat[,k+1]=vMat[,k];
    vMat[s,k+1] = vMat[s,k] + alpha_c * delta;
    #---update state values---
    
    #---update mean threhsolds---
    mMat[,k+1]=mMat[,k];
    mMat[s,k+1] = mMat[s,k] - 
      alpha_m * delta * (mMat[s,k]-positionVect[k]);
    
    if (mMat[s,k+1]<0){
      mMat[s,k+1]=0;
    }
    if (mMat[s,k+1]>300){
      mMat[s,k+1]=300;
    }
    #---update mean threhsolds---
  }
  #===simulate model===
  
  #===compute the likelihood===
  
  if (whatToDo==1){#fitting
    SSE1= sum( ((positionVect[ind1]-mMat[1,ind1])^2)/(2*sigma1^2) );
    SSE2= sum( ((positionVect[ind2]-mMat[2,ind2])^2)/(2*sigma2^2) );
    
    NLL1=sum(log(sigma1))  + SSE1;
    NLL2=sum(log(sigma2))  + SSE2;
    NLL = NLL1+NLL2+.5*numOfTrials*log(2*pi);
  }else if (whatToDo==2){#simulate fitted
    output = list(m1Vect=mMat[1,],m2Vect=mMat[2,])
  }
  #===compute the likelihood===
  
}
