NLL_Model1 = function(x,parameter)
{
  #===parameters===
  #---free parameters---
  if (parameter$numOfThreshold==2){
    sigma1=x[1];sigma2=x[2];
    a1=x[3];a2=x[4];ap1=x[5];ap2=x[6];
    lambda1=x[7];lambda2=x[8];k1=x[9];k2=x[10];
  }else{
    sigma1=x[1];sigma2=sigma1;
    a1=x[2];a2=a1;ap1=x[3];ap2=ap1;
    lambda1=x[4];lambda2=lambda1;k1=x[5];k2=k1;
  }
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
  
  #===simulate model===
  bMat=matrix(0,nrow=2,ncol=numOfTrials);#thrshold in each trial
  bMat[1,ind1]=
    a1-(1-exp(-(rtVect[ind1]/lambda1)^k1))*(.5*a1-ap1);
  bMat[2,ind2]=
    a2-(1-exp(-(rtVect[ind2]/lambda2)^k2))*(.5*a2-ap2);
  #===simulate model===
  
  #===compute the likelihood===
  if (whatToDo==1){#fitting
    n1=sum(ind1);n2=sum(ind2);
    SSE1= sum((positionVect[ind1]-bMat[1,ind1])^2);
    SSE2= sum((positionVect[ind2]-bMat[2,ind2])^2);
    
    NLL1=n1*log(sigma1)  + SSE1/(2*sigma1^2);
    NLL2=n2*log(sigma2)  + SSE2/(2*sigma2^2);
    NLL = NLL1+NLL2+.5*numOfTrials*log(2*pi);
  }else if (whatToDo==2){#simulate fitted
    output = list(b1Vect=bMat[1,],b2Vect=bMat[2,])
  }
  #===compute the likelihood===
  
}
