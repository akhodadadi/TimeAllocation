NLL_Model2 = function(x,parameter){
  
  #===parameters===
  #---free parameters---
  if (parameter$numOfThreshold==2){
    sigma01=x[1];sigma02=x[2];m01=x[3];m02=x[4];
    beta1=x[5];beta2=x[6];alpha1=x[7];alpha2=x[8];
    
    alpha_s1=x[9];alpha_s2=x[10];
  }else{
    sigma1=x[1];sigma2=sigma1;m01=x[2];m02=m01;
    beta1=x[3];beta2=beta1;alpha1=x[4];alpha2=alpha1;
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
  #---compute utility of reward---
  ind=rewardVect>0
  rewardVect[ind]=alpha1*(1-exp(-beta1*rewardVect[ind]))/beta1
  ind=rewardVect<0
  rewardVect[ind]=alpha2*(1-exp(-beta2*rewardVect[ind]))/beta2
  #---compute utility of reward---
  
  #---update---
  m1Vect = matrix(0,nrow=numOfTrials,ncol=1);m1Vect[1]=m01;
  delta1=cumsum(rewardVect*ind1);
  m1Vect[-1] = m01 - delta1[-numOfTrials];
  m1Vect[m1Vect<0]=0;
  m1Vect[m1Vect>300]=300;
  
  m2Vect = matrix(0,nrow=numOfTrials,ncol=1);m2Vect[1]=m02;
  delta2=cumsum(rewardVect*ind2);
  m2Vect[-1] = m02 - delta2[-numOfTrials];
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
