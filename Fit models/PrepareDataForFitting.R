#===determine the directory and load data===
fileName = paste0(dataDir,'SUB',sub,'.csv');

dataMatrix=read.csv(fileName,header=F);
numOfTrials = nrow(dataMatrix);
#===determine the directory and load data===

#===Extract desired columns===
if (sum(c(1,2,3)==ie)!=0){#for experiments A-C
  conditionVect = dataMatrix[,4];
  ind1=conditionVect==1;
  ind2=conditionVect==2;
  rewardVect = dataMatrix[,3];
  rtVect = dataMatrix[,2]+ 4.05;#4.05 sec spent on cue+RSI+...
  positionVect = abs(dataMatrix[,6]);
  
  DPvect=c(3,0);
  rtVect[rewardVect==0 & ind1] = rtVect[rewardVect==0 & ind1]+
    DPvect[1];#add delay penalty to incorrect RT
  rtVect[rewardVect==0 & ind2] = rtVect[rewardVect==0 & ind2]+
    DPvect[2];#add delay penalty to incorrect RT
  
  rewardVect[rewardVect==1 & ind1]=20
  rewardVect[rewardVect==0 & ind1]=-20
  rewardVect[rewardVect==1 & ind2]=1
  rewardVect[rewardVect==0 & ind2]=-1
}else{#experiments D-F
  rewardVect = dataMatrix[,2];
  rtVect = dataMatrix[,3]
  positionVect = abs(dataMatrix[,5]);
  conditionVect = dataMatrix[,7];
  ind1=conditionVect==1;
  ind2=conditionVect==2;
}
#===Extract desired columns===


