NLL_DifficultyDetectionModel3 = function(x,parameter){
  #===parameters===
  #---free parameters---
  sigma1=x[1];sigma2=x[2];
  m10=x[3];m20=x[4];
  alpha_c=x[5];alpha_m1=x[6];alpha_m2=x[7];
  c1=x[8];c2=x[9];
  V10=x[10];V20=x[11];
  t_detect_max=x[12];th_D=x[13];
  #---free parameters---
  
  #---other---
  conditionVect=parameter$conditionVect;
  ind1=conditionVect==1;ind2=conditionVect==2;
  ind=ind1 | ind2
  conditionVect=conditionVect[ind];
  rewardVect=parameter$rewardVect[ind];
  positionVect=parameter$positionVect[ind];
  canoePathMat = parameter$canoePathMat[ind,]
  rtVect=parameter$rtVect[ind];
  ind1=ind1[ind];ind2=ind2[ind];
  numOfTrials=length(rewardVect);
  whatToDo=parameter$whatToDo
  #---other---
  #===parameters===
  
  #===simulate model===
  mMat = matrix(0,nrow=2,ncol=numOfTrials);
  mMat[1,1]=m10;mMat[2,1]=m20;
  
  vMat = matrix(0,nrow=2,ncol=numOfTrials);
  vMat[1,1]=V10;vMat[2,1]=V20;
  
  #---determine detected states in all trials---
  detectedState=rep(2,numOfTrials)
  
  for (k in 1:numOfTrials){
    if (t_detect_max<(rtVect[k]-4.05)){
      canoePath=canoePathMat[k,]
      ind_d = which(abs(canoePath)>=th_D)
      if (length(ind_d)==0){
        t_detect=15;
      }else{
        t_detect=ind_d[1]*1/60;
      }
      
      if (t_detect<=t_detect_max){ #detected that the trial is easy
        detectedState[k]=1;
      }
    }
  }
  #---determine detected states in all trials---
  
  for (k in 1:(numOfTrials-1)){
    #---detected states at k and k+1---
    s=detectedState[k];sPrime=detectedState[k+1];
    #---detected states at k and k+1---
    
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
    ind_d1=detectedState==1;ind_d2=detectedState==2;
    n1=sum(ind_d1);n2=sum(ind_d2);
    SSE1= sum((positionVect[ind_d1]-mMat[1,ind_d1])^2);
    SSE2= sum((positionVect[ind_d2]-mMat[2,ind_d2])^2);
    
    NLL1=n1*log(sigma1)  + SSE1/(2*sigma1^2);
    NLL2=n2*log(sigma2)  + SSE2/(2*sigma2^2);
    NLL = NLL1+NLL2+.5*numOfTrials*log(2*pi);
  }else if (whatToDo==2){#simulate fitted
    output = list(m1Vect=mMat[1,],m2Vect=mMat[2,],
                  detectedState=detectedState)
  }
  #===compute the likelihood===
}
