clear
clc
close all

exper_idx=3;numOfThreshold=1;
method=1;%how to approximate log model evidence. 1:AIC, 2:BIC, 3:Laplace
ms={'AIC','BIC','BIC'};ms=ms{method};%for Laplace we need BIC
%(see Eq.6 in Gershman (2016))


experNames={'A','B','C','D','E','F'};
experiment=experNames{exper_idx};
excludeInd={[5,17,23],3,[1,2]};excludeInd=excludeInd{exper_idx};

%===load BIC/AIC===
rootDir='C:/Users/arakhoda/Dropbox/Projects/Psychology/Learning in SMDP tasks/Experiments/Canoe/';
fn = strcat(rootDir,'AnalyseAndFit/Fitted models/',experiment,'/',ms,'_',num2str(numOfThreshold),'threshold.csv');
[gof_mat,~]=xlsread(fn);
gof_mat(excludeInd,:)=[];gof_mat(:,1)=[];
gof_mat=-gof_mat/2;
%===load BIC/AIC===

%[alpha,exp_r,xp,pxp,bor] = bms(gof_mat(:,[2,9]))%compare model 3 and 8
[alpha,exp_r,xp,pxp,bor] = bms(gof_mat)%compare all models

bar(pxp)