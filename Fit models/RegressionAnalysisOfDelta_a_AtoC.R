#This scripts performs tests on homogeneity of variance. Here
#I use the method proposed in WEISBERG, 2014, CH7.

#===some initial settings===
rm(list=ls())
graphics.off()
rootDir="C:/Users/arakhoda/Dropbox/Projects/Psychology/Learning in SMDP tasks/Experiments/Canoe/"
setwd(rootDir)
source('AnalyseAndFit/Fit models/PartitionData.R')

library(DEoptim)
library(car)
#===some initial settings===

#===specify desired experiment and condition===
desiredSubjects=1:50
ifPerformFullAnalysis=1;#perform hierarchical analysis on 15 models
experiment = 'A'
ifAnalysePooled=0;#if 1, perform analysis on the data pooled
#from all subjects
numOfCond=2
#===specify desired experiment and condition===

#===experiment specific parameters===
expDir=paste0(rootDir,experiment,'/')#experiment's dir
dataDir=paste0(expDir,'Data/')
experNames=c('A','C','D','E','F','SimFromModel3','SimFromModel10')
ie = which(experNames==experiment)
n=c(29,22,11,6,8,100,100)[ie];
if (sum(desiredSubjects>n)){
  desiredSubjects=desiredSubjects[-which(desiredSubjects>n)]  
}
excludeInd=list(c(5,17,23),c(1,2),c(3,9,10),c(2),c(4),0,0)[[ie]];
rhoInEachBlock =
  list(rep(1,40),#Experiment A
       rep(1,35),#Experiment C
       c(rep(1,20),rep(2,15)),#Experiment D
       c(rep(2,10),rep(3,5),rep(4,5),rep(3,5),rep(2,5),rep(1,5)),#Experiment E
       c(rep(2,3),rep(4,3),rep(1,3),rep(3,5),rep(1,4),rep(4,3),
         rep(2,3),rep(3,3),rep(1,3),rep(4,3),rep(2,2)),#Experiment F
       rep(1,40),rep(1,40))
rhoInEachBlock=rhoInEachBlock[[ie]]
n_rho_partitions=c(1,1,2,2,2,1,1)[ie]
n_r_partitions=2
n_t_partitions=c(1,1,1,1,1,1,1)[ie]
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
meanFunc=D_a~cond*r*time
formulaList = list(~cond*r+cond*time,#1
                   ~cond*r+time,#2
                   ~cond*r,#3
                   ~cond*time,#4
                   ~time,~r,~1)#5-7

testPathList=list(c(1,2,3,5),c(1,2,3,6),c(1,2,7),
                  c(1,4,5),c(1,4,7))

formulaList = list(~cond*r+cond*time,~cond*r)
testPathList=list(c(1,2))


subjects_pValues=list()
k=1;
pVal_mat=matrix(0,nrow=numOfSub,ncol=length(formulaList))
chi2_mat=matrix(0,nrow=numOfSub,ncol=length(formulaList))
df_mat=matrix(0,nrow=numOfSub,ncol=length(formulaList))

for (sub in unique(dataFrameForRegression$sub)){
  dataOfThisSub=
    dataFrameForRegression[dataFrameForRegression$sub==sub,]
  m=lm(meanFunc,data=dataOfThisSub)
  compare_pVal_mat=matrix(1,nrow=length(testPathList),ncol=1)
  
  for (formulaCount in 1:length(formulaList)){
    ncv = ncvTest(m,formulaList[[formulaCount]])
    pVal_mat[k,formulaCount]=
      pchisq(ncv$ChiSquare,ncv$Df,lower.tail=FALSE)
    chi2_mat[k,formulaCount]=ncv$ChiSquare
    df_mat[k,formulaCount]=ncv$Df
  }
  
  for (pathCount in 1:length(testPathList)){
    thisPath=testPathList[[pathCount]]
    pathLength=length(thisPath)-1
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

#acf(dataFrameForRegression$D_a,lag.max = 5)

#---extract simplest model of variance---
simplestModelOfVar=matrix(0,nrow=numOfSub,
                          ncol=nrow(compare_pVal_mat))
for (sub in 1:numOfSub){
  compare_pVal_mat = subjects_pValues[[sub]]
  ind=compare_pVal_mat<.05
  for (pathCount in 1:nrow(compare_pVal_mat)){
    if (sum(ind[pathCount,])){
      simplestModelOfVar[sub,pathCount]=
        testPathList[[pathCount]][which(ind[pathCount,])[1]]
    }else{
      simplestModelOfVar[sub,pathCount]=
        testPathList[[pathCount]][ncol(compare_pVal_mat)]
    }
  }
  simplestModelOfVar[sub,pVal_mat[sub,simplestModelOfVar[sub,]]>.05]=0
  
  #---non-nested models---
  #compare model3 and model4
  ch2_1=chi2_mat[sub,3]
  ch2_2=chi2_mat[sub,4]
  df1=df_mat[sub,3]
  df2=df_mat[sub,4]
  if (ch2_2<ch2_1){
    ch2=pchisq(ch2_1-ch2_2,df1-df2,lower.tail=FALSE)
    if (ch2<.05){
      simplestModelOfVar[simplestModelOfVar==4]=3
    }
  }
  #compare model3 and model4
  
  #compare model2 and model4
  ch2_1=chi2_mat[sub,2]
  ch2_2=chi2_mat[sub,4]
  df1=df_mat[sub,2]
  df2=df_mat[sub,4]
  if (ch2_2<ch2_1){
    ch2=pchisq(ch2_1-ch2_2,df1-df2,lower.tail=FALSE)
    if (ch2<.05){
      simplestModelOfVar[simplestModelOfVar==4]=2
    }
  }
  #compare model2 and model4
  
  #compare model2 and model3
  ch2_1=chi2_mat[sub,2]
  ch2_2=chi2_mat[sub,3]
  df1=df_mat[sub,2]
  df2=df_mat[sub,3]
  if (ch2_2<ch2_1){
    ch2=pchisq(ch2_1-ch2_2,df1-df2,lower.tail=FALSE)
    if (ch2<.05){
      simplestModelOfVar[simplestModelOfVar==3]=2
    }
  }
  #compare model2 and model3
  #---non-nested models---
}  

m1=apply(simplestModelOfVar,1,max)
windows()
hist(m1,breaks = seq(-.5,7.5,1))
#---extract simplest model of variance---

#---for SimFromModel3, compare models 1 and 5---
if (ie>=6){
  compare1vs5_pVal=rep(0,numOfSub)
  for (sub in 1:numOfSub){
    ch2_1=chi2_mat[sub,2]
    ch2_5=chi2_mat[sub,5]
    df1=df_mat[sub,2]
    df5=df_mat[sub,5]
    compare1vs5_pVal[sub]=
      pchisq(ch2_1-ch2_5,df1-df5,lower.tail=FALSE)
  }
  ind=compare1vs5_pVal<.1
  chosenModel=rep(1,numOfSub)
  chosenModel[ind]=2#model 1
  windows();hist(chosenModel)
}
#---for SimFromModel3, compare models 1 and 5---
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%ANALYSE HOMOGENEITY OF VARIANCE%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

