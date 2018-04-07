NLL_Model5 = function(x,parameter){
  
  #===parameters===
  #---free parameters---
  if (parameter$numOfThreshold==2){
    sigma01=x[1];sigma02=x[2];m01=x[3];m02=x[4];
    a1=x[5];a2=x[6];b1=x[7];b2=x[8];
    theta1=x[9];theta2=x[10];
    alpha1=x[11];alpha2=x[12];
    
    alpha_s1=x[13];alpha_s2=x[14];
  }else{
    sigma1=x[1];sigma2=sigma1;m01=x[2];m02=m01;
    a1=x[3];a2=a1;b1=x[4];b2=b1;
    theta1=x[5];theta2=theta1;
    alpha1=x[6];alpha2=alpha1;
  }
  #---free parameters---
  
  #---other---
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
  
  #---update---
  m1Vect = matrix(0,nrow=numOfTrials,ncol=1);m1Vect[1]=m01;
  delta1=cumsum(pHatVect[ind1]-theta1);delta1=c(0,delta1[-1]);
  m1Vect[ind1] = m01 - alpha1*delta1;
  m1Vect[m1Vect<0]=0;
  m1Vect[m1Vect>300]=300;
  
  m2Vect = matrix(0,nrow=numOfTrials,ncol=1);m2Vect[1]=m02;
  delta2=cumsum(pHatVect[ind2]-theta2);delta2=c(0,delta2[-1]);
  m2Vect[ind2] = m02 - alpha2*delta2;
  m2Vect[m2Vect<0]=0;
  m2Vect[m2Vect>300]=300;
  #---update---
  #===simulate model===
  
  #===compute the likelihood===
  if (whatToDo==1){#fitting
    SSE1= sum( ((positionVect[ind1]-m1Vect[ind1])^2)/(2*sigma1^2) );
    SSE2= sum( ((positionVect[ind2]-m2Vect[ind2])^2)/(2*sigma2^2) );
    
    NLL1=sum(log(sigma1))  + SSE1;
    NLL2=sum(log(sigma2))  + SSE2;
    NLL = NLL1+NLL2+.5*numOfTrials*log(2*pi);
  }else if (whatToDo==2){#simulate fitted
    output = list(m1Vect=m1Vect,m2Vect=m2Vect)
  }
  #===compute the likelihood===
}
