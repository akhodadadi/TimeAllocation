#This scripts performs tests on homogeneity of variance. Here
#I use the method proposed in WEISBERG, 2014, CH7.

#===some initial settings===
rm(list=ls())
graphics.off()
rootDir="C:/Users/arakhoda/Dropbox/Projects/Psychology/Learning in SMDP tasks/Experiments/Canoe/"
setwd(rootDir)
source('AnalyseAndFit/Models/SourceAllModels_UsingDelta_a.R')
source('AnalyseAndFit/Fit models/PartitionData.R')

library(DEoptim)
library(car)
#===some initial settings===

#===specify desired experiment and condition===
desiredSubjects=1
ifPerformFullAnalysis=1;#perform hierarchical analysis on 15 models
experiment = 'A'
ifAnalysePooled=0;#if 1, perform analysis on the data pooled
#from all subjects
#===specify desired experiment and condition===

#===experiment specific parameters===
expDir=paste0(rootDir,experiment,'/')#experiment's dir
dataDir=paste0(expDir,'Data/')
experNames=c('A','C','D','E','F')
ie = which(experNames==experiment)
n=c(29,22,11,6,8)[ie];
excludeInd=list(c(5,17,23),c(1,2),c(3,9,10),c(2),c(4))[[ie]];
numOfParams = numOfFreeParams[[modelNo]]
rhoInEachBlock =
  list(rep(1,40),#Experiment A
       rep(1,35),#Experiment C
       c(rep(1,20),rep(2,15)),#Experiment D
       c(rep(2,10),rep(3,5),rep(4,5),rep(3,5),rep(2,5),rep(1,5)),#Experiment E
       c(rep(2,3),rep(4,3),rep(1,3),rep(3,5),rep(1,4),rep(4,3),
         rep(2,3),rep(3,3),rep(1,3),rep(4,3),rep(2,2)) )#Experiment F
rhoInEachBlock=rhoInEachBlock[[ie]]
n_rho_partitions=c(1,1,2,2,2)[ie]
n_r_partitions=2
n_t_partitions=c(1,1,1,1,1)[ie]
#===experiment specific parameters===

#===prepare data===
dataFrameForRegression = data.frame()
for (sub in desiredSubjects){
  cat("\014")
  print(sub)
  if (sum(excludeInd==sub))
    next;
  
  source("AnalyseAndFit/Fit models/PrepareDataForFitting_usingDelta_a.R")
  d=data.frame(rep(sub,nrow(partitionedDataFrame)));
  colnames(d)="sub"
  d=cbind(d,partitionedDataFrame)
  dataFrameForRegression=
    rbind(dataFrameForRegression,d)
}
if (ifAnalysePooled){
  dataFrameForRegression$sub=1
}
numOfSub=length(unique(dataFrameForRegression$sub))
#===prepare data===


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%ANALYSE HOMOGENEITY OF VARIANCE%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (ie<=2){#experiments A,C
  meanFunc=D_a~cond*r*time
  formulaList = list(~cond*r*time,#1
                     ~cond+r+time+cond*r+r*time+cond*time,#2
                     ~cond+r+time+cond*r+r*time,#3
                     ~cond+r+time+cond*r+cond*time,#4
                     ~cond+r+time+r*time+cond*time,#5
                     ~cond+r+time+cond*r,#6
                     ~cond+r+time+r*time,#7
                     ~cond+r+time+cond*time,#8
                     ~cond+r+time,#9
                     ~cond+r,~cond+time,~r+time,#10-12
                     ~cond,~r,~time)#13-15
}
if (ie>2){#experiments D-F
  meanFunc=D_a~cond*r*time*rho
  formulaList = list(~cond*r*time*rho,#1
                     ~cond+r+time*rho+cond*r+time*r*rho+time*cond*rho,#2
                     ~cond+r+time*rho+cond*r+r*time*rho,#3
                     ~cond+r+time*rho+cond*r+cond*time*rho,#4
                     ~cond+r+time*rho+r*time*rho+cond*time*rho,#5
                     ~cond+r+time*rho+cond*r,#6
                     ~cond+r+time*rho+r*time*rho,#7
                     ~cond+r+time*rho+cond*time*rho,#8
                     ~cond+r+time*rho,#9
                     ~cond+time*r,~cond+time*rho,~r+time*rho,#10-12
                     ~cond,~r,~time*rho)#13-15
}

#%%%%%%full analysis%%%%%%%%%
if (ifPerformFullAnalysis){
  testPathList=list(c(1,2,3,6,9,10,13),c(1,2,3,6,9,10,14),
                    c(1,2,3,6,9,11,13),c(1,2,3,6,9,11,15),
                    c(1,2,3,6,9,12,14),c(1,2,3,6,9,12,15),
                    
                    c(1,2,4,6,9,10,13),c(1,2,4,6,9,10,14),
                    c(1,2,4,6,9,11,13),c(1,2,4,6,9,11,15),
                    c(1,2,4,6,9,12,14),c(1,2,4,6,9,12,15),
                    
                    c(1,2,3,7,9,10,13),c(1,2,3,7,9,10,14),
                    c(1,2,3,7,9,11,13),c(1,2,3,7,9,11,15),
                    c(1,2,3,7,9,12,14),c(1,2,3,7,9,12,15),
                    
                    c(1,2,5,7,9,10,13),c(1,2,5,7,9,10,14),
                    c(1,2,5,7,9,11,13),c(1,2,5,7,9,11,15),
                    c(1,2,5,7,9,12,14),c(1,2,5,7,9,12,15),
                    
                    c(1,2,4,8,9,10,13),c(1,2,4,8,9,10,14),
                    c(1,2,4,8,9,11,13),c(1,2,4,8,9,11,15),
                    c(1,2,4,8,9,12,14),c(1,2,4,8,9,12,15),
                    
                    c(1,2,5,8,9,10,13),c(1,2,5,8,9,10,14),
                    c(1,2,5,8,9,11,13),c(1,2,5,8,9,11,15),
                    c(1,2,5,8,9,12,14),c(1,2,5,8,9,12,15))
  
  pathLength=6;
  subjects_pValues=list()
  k=1;
  pVal_mat=matrix(0,nrow=numOfSub,ncol=length(formulaList))
  chi2_mat=matrix(0,nrow=numOfSub,ncol=length(formulaList))
  df_mat=matrix(0,nrow=numOfSub,ncol=length(formulaList))
  
  for (sub in unique(dataFrameForRegression$sub)){
    dataOfThisSub=
      dataFrameForRegression[dataFrameForRegression$sub==sub,]
    m=lm(meanFunc,data=dataOfThisSub)
    compare_pVal_mat=matrix(0,nrow=36,ncol=pathLength)
    
    for (formulaCount in 1:length(formulaList)){
      ncv = ncvTest(m,formulaList[[formulaCount]])
      pVal_mat[k,formulaCount]=
        pchisq(ncv$ChiSquare,ncv$Df,lower.tail=FALSE)
      chi2_mat[k,formulaCount]=ncv$ChiSquare
      df_mat[k,formulaCount]=ncv$Df
    }
    
    for (pathCount in 1:length(testPathList)){
      thisPath=testPathList[[pathCount]]
      for (testCount in 1:pathLength){
        ch2_1=chi2_mat[k,thisPath[testCount]]
        ch2_2=chi2_mat[k,thisPath[testCount+1]]
        df1=df_mat[k,thisPath[testCount]]
        df2=df_mat[k,thisPath[testCount+1]]
        compare_pVal_mat[pathCount,testCount]=
          pchisq(ch2_1-ch2_2,df1-df2,lower.tail=FALSE)
      }
    }
    subjects_pValues[[k]]=compare_pVal_mat
    k=k+1
  }
  View(compare_pVal_mat)
  View(pVal_mat)
}
#%%%%%%full analysis%%%%%%%%%


#%%%%%%short analysis%%%%%%%%%
pVal_short=matrix(0,nrow=numOfSub,ncol=6)
k=1;
for (sub in unique(dataFrameForRegression$sub)){
  dataOfThisSub=
    dataFrameForRegression[dataFrameForRegression$sub==sub,]
  m=lm(meanFunc,data=dataOfThisSub)
  if (ie<=2){#experiments A,C
    ncv_full=ncvTest(m,~r*cond*time)
    ncv_model3=ncvTest(m,~r*cond+time)
    ncv_model10=ncvTest(m,~cond)
  }
  if (ie>2){#experiments D-F
    ncv_full=ncvTest(m,~r*cond*time*rho)
    ncv_model3=ncvTest(m,~r*cond+time*rho)
    ncv_model10=ncvTest(m,~cond)
  }
  pVal_short[k,1]=
    pchisq(ncv_full$ChiSquare,ncv_full$Df,lower.tail=FALSE)
  pVal_short[k,2]=
    pchisq(ncv_model3$ChiSquare,ncv_model3$Df,lower.tail=FALSE)
  pVal_short[k,3]=
    pchisq(ncv_model10$ChiSquare,ncv_model10$Df,lower.tail=FALSE)
  pVal_short[k,4]=
    pchisq(ncv_full$ChiSquare-ncv_model3$ChiSquare,
           ncv_full$Df-ncv_model3$Df,lower.tail=FALSE)
  pVal_short[k,5]=
    pchisq(ncv_full$ChiSquare-ncv_model10$ChiSquare,
           ncv_full$Df-ncv_model10$Df,lower.tail=FALSE)
  pVal_short[k,6]=
    pchisq(ncv_model3$ChiSquare-ncv_model10$ChiSquare,
           ncv_model3$Df-ncv_model10$Df,lower.tail=FALSE)
  k=k+1
}

View(pVal_short)
#%%%%%%short analysis%%%%%%%%%
#===regression analysis===

acf(dataFrameForRegression$D_a,lag.max = 5)

#===extract simplest model of variance===
if (ifPerformFullAnalysis==1){
  simplestModelOfVar=matrix(0,nrow=numOfSub,
                            ncol=nrow(compare_pVal_mat))
  for (sub in 1:numOfSub){
    compare_pVal_mat = subjects_pValues[[sub]]
    ind=compare_pVal_mat<.05
    for (pathCount in 1:nrow(compare_pVal_mat)){
      if (sum(ind[pathCount,])){
        simplestModelOfVar[sub,pathCount]=
          testPathList[[pathCount]][which(ind[pathCount,])[1]]
      }
    }
    simplestModelOfVar[sub,pVal_mat[sub,simplestModelOfVar[sub,]]>.05]=0
  }  
}

apply(simplestModelOfVar,1,sum)
apply(pVal_mat<.05,1,sum)
m1=apply(simplestModelOfVar,1,max)
hist(m1,breaks = seq(-.5,15.5,1))
#===extract simplest model of variance===
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%ANALYSE HOMOGENEITY OF VARIANCE%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%ANALYSE MEAN%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%ANALYSE MEAN%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


