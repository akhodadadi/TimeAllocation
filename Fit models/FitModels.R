cat('\014')
rm(list=ls())
graphics.off()
rootDir="C:/Users/arakhoda/Dropbox/Projects/Psychology/Learning in SMDP tasks/Experiments/Canoe/"
setwd(rootDir)
library(DEoptim)

#===specify desired experiment/model/condition===
modelNo=11
experiment = 'C'
numOfThreshold=1

source('AnalyseAndFit/Models/SourceAllModels.R')
#===specify desired experiment/model/condition===

expDir=paste0(rootDir,experiment,'/')#experiment's dir
dataDir=paste0(expDir,'Data/')
experNames=c('A','B','C','D','E','F','SimFromModel3','SimFromModel10')
ie = which(experNames==experiment)
n=c(29,2,22,11,6,8,50,50)[ie];
excludeInd=list(c(5,17,23),0,c(1,2),c(3,9,10),c(2),c(4),0,0)[[ie]];

numOfParams = numOfFreeParams[[modelNo+1]]
fittedParamsMat = matrix(0,nrow = n,ncol = numOfParams)
NLL = matrix(0,nrow = n,ncol = 1)
fittedModelsList = list();
N = matrix(0,nrow = n,ncol = 1);#no. of trials

for (sub in 1:n){
  cat("\014")
  print(sub)
  if (sum(excludeInd==sub))
    next;
  
  #===load data===
  source("AnalyseAndFit/Fit models/PrepareDataForFitting.R")
  #===load data===
  
  #===fit the model===
  if (modelNo==0){
    if (numOfThreshold==2){
      sigma1=sd(positionVect[ind1]);m1=mean(positionVect[ind1])
      sigma2=sd(positionVect[ind2]);m2=mean(positionVect[ind2])
      optim=list(bestmem=c(sigma1,sigma2,m1,m2),bestval=0)
    }else{
      sigma1=sd(positionVect[ind1 | ind2]);
      m1=mean(positionVect[ind1 | ind2])
      sigma2=sigma1;m2=m1;
      optim=list(bestmem=c(sigma1,m1),bestval=0)
    }
    SSE1=sum((positionVect[ind1]-m1)^2);n1=sum(ind1)
    SSE2=sum((positionVect[ind2]-m2)^2);n2=sum(ind2)
    NLL1 = (n1*log(sigma1)) + (SSE1/(2*sigma1^2));
    NLL2 = (n2*log(sigma2)) + (SSE2/(2*sigma2^2));
    optim$bestval=NLL1+NLL2+.5*(n1+n2)*log(2*pi);
    D=list(optim=optim)

  }else{
    parameter = list(rewardVect=rewardVect,
                     positionVect=positionVect,
                     rtVect=rtVect,
                     conditionVect=conditionVect,
                     whatToDo=1,
                     numOfThreshold=numOfThreshold);
    D = DEoptim(modelsName[[modelNo+1]],
                lower=lowerLimit[[modelNo+1]],
                upper=upperLimit[[modelNo+1]],
                control=list(NP=10*numOfParams+30,
                             itermax=500,trace=1),
                parameter=parameter)
  }
  #===fit the model===
  
  #===save fitted model for each subject===
  fittedParamsMat[sub,] = D$optim$bestmem
  NLL[sub,] = D$optim$bestval
  fittedModelsList[[sub]]= D
  N[sub,] = sum(ind1 | ind2)
  #===save fitted model for each subject===
}

#===save fitted parameters in csv file===
fn = paste0('AnalyseAndFit/Fitted models/',
            experiment,'/model',modelNo,'_',
            numOfThreshold,'threshold_4','.csv')
d=data.frame(fittedParamsMat,NLL,N);
colnames(d)=col.names=c(paramNames[[modelNo+1]],'NLL','N')
write.csv(d,file=fn)
#===save fitted parameters in csv file===

#===save fitted parameters in Rdata file===
# fn = paste0('AnalyseAndFit/Fitted models/',
#             experiment,'/model',modelNo,'_',
#             numOfThreshold,'threshold_1','.Rdata')
# save(fittedModelsList,file=fn)
#===save fitted parameters in Rdata file===


