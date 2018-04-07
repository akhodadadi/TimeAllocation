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
desiredSubjects=1:11
ifPerformFullAnalysis=1;#perform hierarchical analysis on 15 models
experiment = 'F'
ifAnalysePooled=0;#if 1, perform analysis on the data pooled
#from all subjects
numOfCond=2
#===specify desired experiment and condition===

#===experiment specific parameters===
expDir=paste0(rootDir,experiment,'/')#experiment's dir
dataDir=paste0(expDir,'Data/')
experNames=c('A','C','D','E','F')
ie = which(experNames==experiment)
n=c(29,22,11,6,8)[ie];
if (sum(desiredSubjects>n)){
  desiredSubjects=desiredSubjects[-which(desiredSubjects>n)]  
}
excludeInd=list(c(5,17,23),c(1,2),c(3,9,10),c(2),c(4))[[ie]];
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
if (numOfCond==1){n_r_partitions=2*n_r_partitions}
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
meanFunc=c(D_a~r*time*rho,D_a~cond*r*time*rho)[[numOfCond]]
# formulaList = list(~cond*r+cond*time*rho,#1
#                    ~cond*r+time*rho,#2
#                    ~cond*r+time*cond,#3
#                    ~cond*r+time,#4
#                    ~cond,~r,~time)#5-7

formulaList = list(~cond*r+time*rho,#1
                   ~cond*r+time,#2
                   ~time*cond)#3
if (numOfCond==1){
  formulaList = list(~r+time*rho,#1
                     ~r+time,#2
                     ~r)#3
}

nModels=length(formulaList)

subjects_pValues=list()
k=1;
pVal_mat=matrix(0,nrow=numOfSub,ncol=length(formulaList))
chi2_mat=matrix(0,nrow=numOfSub,ncol=length(formulaList))
df_mat=matrix(0,nrow=numOfSub,ncol=length(formulaList))
compare_pVal_mat=matrix(1,nrow=numOfSub,ncol=nModels)

for (sub in unique(dataFrameForRegression$sub)){
  dataOfThisSub=
    dataFrameForRegression[dataFrameForRegression$sub==sub,]
  m=lm(meanFunc,data=dataOfThisSub)
  
  for (formulaCount in 1:nModels){
    ncv = ncvTest(m,formulaList[[formulaCount]])
    pVal_mat[k,formulaCount]=
      pchisq(ncv$ChiSquare,ncv$Df,lower.tail=FALSE)
    chi2_mat[k,formulaCount]=ncv$ChiSquare
    df_mat[k,formulaCount]=ncv$Df
  }
  
  ch2_1=chi2_mat[k,1];df1=df_mat[k,1]
  ch2_2=chi2_mat[k,2];df2=df_mat[k,2]
  ch2_3=chi2_mat[k,3];df3=df_mat[k,3]
  compare_pVal_mat[k,1]=pchisq(abs(ch2_1-ch2_2),abs(df1-df2),
                               lower.tail=FALSE)
  compare_pVal_mat[k,2]=pchisq(abs(ch2_1-ch2_3),abs(df1-df3),
                               lower.tail=FALSE)
  compare_pVal_mat[k,3]=pchisq(abs(ch2_2-ch2_3),abs(df2-df3),
                               lower.tail=FALSE)
  k=k+1
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%ANALYSE HOMOGENEITY OF VARIANCE%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


View(compare_pVal_mat)
View(pVal_mat)


