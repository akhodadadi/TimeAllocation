NLL_Model6 = function(x,parameter){
  
  #===parameters===
  #---free parameters---
  if (parameter$numOfThreshold==2){
    sigma01=x[1];sigma02=x[2];
    a01=x[3];a02=x[4];ap1=x[5];ap2=x[6];
    lambda1=x[7];lambda2=x[8];k1=x[9];k2=x[10];
    
    a1=x[11];a2=x[12];b1=x[13];b2=x[14];
    theta1=x[15];theta2=x[16];
    alpha1=x[17];alpha2=x[18];
    
    alpha_s1=x[19];alpha_s2=x[20];
  }else{
    
  }
  #---free parameters---
  
  #---other---
  rtVect=parameter$rtVect;
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
  
  #===simulate model===
  #---estimate accuracy---
  pHatVect = matrix(0,nrow=numOfTrials,ncol=1);
  trials1=1:sum(ind1);trials2=1:sum(ind2);
  ifCorrect1=rewardVect[ind1]>0;ifCorrect2=rewardVect[ind2]>0;
  pHatVect[ind1]=(a1+cumsum(ifCorrect1))/(a1+b1+trials1);
  pHatVect[ind2]=(a2+cumsum(ifCorrect2))/(a2+b2+trials2);
  #---estimate accuracy---
  
  #---update 'a' in Weibull---
  m1Vect = matrix(0,nrow=numOfTrials,ncol=1);
  delta1=cumsum(pHatVect[ind1]-theta1);delta1=c(0,delta1[-1]);
  m1Vect[ind1] = -alpha1*delta1;
  m1Vect[m1Vect+a01<0]=-a01;
  m1Vect[m1Vect+a01>300]=300-a01;
  
  m2Vect = matrix(0,nrow=numOfTrials,ncol=1);
  delta2=cumsum(pHatVect[ind2]-theta2);delta2=c(0,delta2[-1]);
  m2Vect[ind2] = -alpha2*delta2;
  m2Vect[m2Vect+a02<0]=-a02;
  m2Vect[m2Vect+a02>300]=300-a02;
  #---update 'a' in Weibull---
  
  #---compute Weibull threshold in each trial---
  a1=m1Vect[ind1];a2=m2Vect[ind2];
  bMat=matrix(0,nrow=2,ncol=numOfTrials);#thrshold in each trial
  bMat[1,ind1]=
    a1+a01-(1-exp(-(rtVect[ind1]/lambda1)^k1))*(.5*a01-ap1);
  bMat[2,ind2]=
    a2+a02-(1-exp(-(rtVect[ind2]/lambda2)^k2))*(.5*a02-ap2);
  #---compute Weibull threshold in each trial---
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
