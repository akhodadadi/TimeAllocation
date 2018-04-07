NLL_Model7 = function(x,parameter){
  
  #===parameters===
  #---free parameters---
  if (parameter$numOfThreshold==2){
    sigma01=x[1];sigma02=x[2];
    a01=x[3];a02=x[4];ap1=x[5];ap2=x[6];
    lambda1=x[7];lambda2=x[8];k1=x[9];k2=x[10];
    
    r_win=floor(x[11]);t_win=floor(x[12]);
    alpha1=x[13];alpha2=x[14];
    c1=x[15];c2=x[16];
    
    alpha_s1=x[17];alpha_s2=x[18];
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
  
  #===rt and reward for each condition===
  correctVect=rep(0,numOfTrials);correctVect[rewardVect>0]=1;
  correctVect1=correctVect[ind1];
  correctVect2=correctVect[ind2];
  rtVect1=rtVect[ind1];rtVect2=rtVect[ind2];
  #===rt and reward for each condition===
  
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
  
  i1=0;i2=0;
  
  for (k in 1:(numOfTrials-1)){
    s=conditionVect[k]
    
    #---set parameters for current trial---
    if (s==1){
      i1=i1+1;rt=rtVect1;corrVect=correctVect1;i0=i1;c=c1;
      a=aMat[1,k];ap=ap1;lambda=lambda1;k_b=k1;#threshold params
    }
    if (s==2){
      i2=i2+1;rt=rtVect2;corrVect=correctVect2;i0=i2;c=c2;
      a=aMat[2,k];ap=ap2;lambda=lambda2;k_b=k2;#threshold params
    }
    #---set parameters for current trial---
    
    #---estimate accuracy---
    if (i0>r_win){
      p_hat = sum(corrVect[(i0-r_win):i0])/r_win
    }else if (i0<=r_win){
      p_hat = sum(corrVect[1:i0])/i0
    }
    #---estimate accuracy---
    
    #---estimate mean RT---
    if (i0>t_win){
      t_hat = sum(rt[(i0-t_win):i0])/t_win
    }else if (i0<=t_win){
      t_hat = sum(rt[1:i0])/i0
    }
    #---estimate mean RT---
    
    #---compute delta---
    delta=alpha1*p_hat+alpha2*t_hat+c
    #---compute delta---
    
    #---compute value of threshold---
    bMat[s,k] = a-(1-exp(-(rtVect[k]/lambda)^k_b))*(.5*a-ap);
    #---compute value of threshold---
    
    #---update 'a' in Weibull threhsolds---
    aMat[,k+1]=aMat[,k];
    aMat[s,k+1] = aMat[s,k] + delta
    
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