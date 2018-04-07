cat('\014')
rm(list=ls())
graphics.off()
rootDir="C:/Users/arakhoda/Dropbox/Projects/Psychology/Learning in SMDP tasks/Experiments/Canoe/"
setwd(rootDir)
library(DEoptim)

#===specify desired experiment/model/condition===
experiment = 'C'
modelVerstion=3;#1:constant cost of time.2:rho as cost of time
                #3:t_D is a gamma R.V.
source('AnalyseAndFit/Models/NLL_DifficultyDetectionModel1.R')
source('AnalyseAndFit/Models/NLL_DifficultyDetectionModel2.R')
source('AnalyseAndFit/Models/NLL_DifficultyDetectionModel3.R')

modelName=list(NLL_DifficultyDetectionModel1,
               NLL_DifficultyDetectionModel2,
               NLL_DifficultyDetectionModel3)[[modelVerstion]]
expDir=paste0(rootDir,experiment,'/')#experiment's dir
dataDir=paste0(expDir,'Data/')
experNames=c('A','C','D','E','F')
ie = which(experNames==experiment)
n=c(29,22,11,6,8,50,50)[ie];
excludeInd=list(c(5,17,23),c(1,2),c(3,9,10),c(2),c(4),0,0)[[ie]];
#===specify desired experiment/model/condition===

#===determine model and experiment dependent parameters===
numOfParams = c(13,13,13)[modelVerstion]
fittedParamsMat = matrix(0,nrow = n,ncol = numOfParams)
NLL = matrix(0,nrow = n,ncol = 1)
AIC = matrix(0,nrow = n,ncol = 1)
BIC = matrix(0,nrow = n,ncol = 1)
fittedModelsList = list();
N = matrix(0,nrow = n,ncol = 1);#no. of trials
paramNames=list(
  c('sigma1','sigma2','m0_1','m0_2','alpha_c','alpha_m1',
    'alpha_m2','c1','c2','v0_1','v0_2','t_detect_max','th_D'),
  c('sigma1','sigma2','m0_1','m0_2','alpha_c','alpha_m1',
    'alpha_m2','rhowin','v0_1','v0_2','v0_3','t_detect_max','th_D'),
  c('sigma1','sigma2','m0_1','m0_2','alpha_c','alpha_m1',
    'alpha_m2','c1','c2','v0_1','v0_2','t_detect_max','th_D'))[[modelVerstion]]
lowerLimit=list(c(0,0,0,0,-.01,-.1,-.1,-20,-20,-15,-15,0,0),
                c(0,0,0,0,-.01,-.1,-.1,1,-15,-15,-15,0,0),
                c(0,0,0,0,-.01,-.1,-.1,-20,-20,-15,-15,0,0))[[modelVerstion]]
upperLimit=list(c(100,100,300,300,.01,.1,.1,20,20,15,15,12,300),
                c(100,100,300,300,.01,.1,.1,50,15,15,15,12,300),
                c(100,100,300,300,.01,.1,.1,20,20,15,15,12,300))[[modelVerstion]]
#===determine model and experiment dependent parameters===

#===model fitting===
for (sub in 17:22){
  cat("\014")
  print(sub)
  if (sum(excludeInd==sub))
    next;
  
  #===load data===
  source("AnalyseAndFit/Fit models/PrepareDataForFitting.R")
  
  #---canoe path---
  fn = paste0(dataDir,'sub',sub,'_CanoeFullPath.csv')
  canoePathMat=scan(fn,what=list(pos=0),sep=",");
  canoePathMat = matrix(canoePathMat$pos,
                        nrow=numOfTrials,
                        ncol=length(canoePathMat$pos)/numOfTrials,
                        byrow=T)
  #---canoe path---
  #===load data===
  
  #===estimate parameters===
  parameter = list(rewardVect=rewardVect,
                   positionVect=positionVect,
                   rtVect=rtVect,
                   conditionVect=conditionVect,
                   whatToDo=1,
                   canoePathMat=canoePathMat);
  
  D = DEoptim(modelName,
              lower=lowerLimit,
              upper=upperLimit,
              control=list(NP=150,itermax=500,trace=1),
              parameter=parameter)
  #===estimate parameters===
  
  #===save fitted model for each subject===
  N[sub,] = sum(ind1 | ind2)
  NLL[sub,] = D$optim$bestval
  fittedParamsMat[sub,] = D$optim$bestmem
  fittedModelsList[[sub]]= D
  
  #---compute goodness of fit---
  AIC[sub,]=2*NLL[sub,] + 2*numOfParams
  BIC[sub,]=2*NLL[sub,] + log(N[sub,])*numOfParams
  #---compute goodness of fit---
  #===save fitted model for each subject===
  
  #===save fitted parameters in csv file===
  fn = paste0('AnalyseAndFit/Fitted models/',
              experiment,'/DiffDetectModel_',
              modelVerstion,'_4.csv')
  d=data.frame(fittedParamsMat,N,NLL,AIC,BIC);
  colnames(d)=col.names=c(paramNames,
                          'N','NLL','AIC','BIC')
  write.csv(d,file=fn)
  #===save fitted parameters in csv file===
}
#===model fitting===



#===save fitted parameters in Rdata file===
# fn = paste0('AnalyseAndFit/Fitted models/',
#             experiment,'/DiffDetectModel_',
#             modelVerstion,'_1.Rdata')
# save(fittedModelsList,file=fn)
#===save fitted parameters in Rdata file===








