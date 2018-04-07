NLL_Model3 = function(x,parameter){
  
  #===parameters===
  #---free parameters---
  if (parameter$numOfThreshold==2){
    sigma01=x[1];sigma02=x[2];
    a01=x[3];a02=x[4];ap1=x[5];ap2=x[6];
    lambda1=x[7];lambda2=x[8];k1=x[9];k2=x[10];
    beta1=x[11];beta2=x[12];alpha1=x[13];alpha2=x[14];
    
    alpha_s1=x[15];alpha_s2=x[16];
  }else{

  }
  #---free parameters---
  
  #---other---
  conditionVect=parameter$conditionVect;
  rewardVect=parameter$rewardVect;
  positionVect=parameter$positionVect;
  numOfTrials=length(rewardVect);
  rtVect=parameter$rtVect;
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
  
  #---update 'a' in Weibull---
  m1Vect = matrix(0,nrow=numOfTrials,ncol=1);
  delta1=cumsum(rewardVect*ind1);
  m1Vect[-1] = -delta1[-numOfTrials];
  m1Vect[m1Vect+a01<0]=-a01;
  m1Vect[m1Vect+a01>300]=300-a01;
  
  m2Vect = matrix(0,nrow=numOfTrials,ncol=1);
  delta2=cumsum(rewardVect*ind2);
  m2Vect[-1] = -delta2[-numOfTrials];
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
  }else if (whatToDo==3){#output param. of threshold
    threshold_param_mat=cbind(a01+m1Vect,
                              rep(ap1,numOfTrials),
                              rep(lambda1,numOfTrials),
                              rep(k1,numOfTrials),
                              a02+m2Vect,
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
