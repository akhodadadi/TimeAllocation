source('AnalyseAndFit/Models/NLL_Model1.r')
source('AnalyseAndFit/Models/NLL_Model2.r')
source('AnalyseAndFit/Models/NLL_Model3.r')
source('AnalyseAndFit/Models/NLL_Model4.r')
source('AnalyseAndFit/Models/NLL_Model5.r')
source('AnalyseAndFit/Models/NLL_Model6.r')
source('AnalyseAndFit/Models/NLL_Model7.r')
source('AnalyseAndFit/Models/NLL_Model8.r')
source('AnalyseAndFit/Models/NLL_Model9.r')
source('AnalyseAndFit/Models/NLL_Model10.r')
source('AnalyseAndFit/Models/NLL_Model11.r')
source('AnalyseAndFit/Models/NLL_Model12.r')

modelsName = c('model0',NLL_Model1,NLL_Model2,NLL_Model3,
               NLL_Model4,NLL_Model5,NLL_Model6,NLL_Model7,
               NLL_Model8,NLL_Model9,NLL_Model10,NLL_Model11,
               NLL_Model12)

if (numOfThreshold==2){
  numOfFreeParams = c(4,10,10,16,18,14,20,18,13,19,20,0,0)
  paramNames = list(c('sigma1','sigma2','m01','m02'),#model0
                    c('sigma1','sigma2','a01','a02','ap1','ap2','lambda1','lambda2','k1','k2'),#model1
                    c('sigma1','sigma2','m01','m02','beta1','beta2','alpha1','alpha2','alpha_s1','alpha_s2'),#model2
                    c('sigma1','sigma2','a01','a02','ap1','ap2','lambda1','lambda2','k1','k2','beta1','beta2','alpha1','alpha2','alpha_s1','alpha_s2'),#model3
                    c('sigma1','sigma2','a01','a02','ap1','ap2','lambda1','lambda2','k1','k2','beta1','beta2','alpha1','alpha2','gamma1','gamma2','alpha_s1','alpha_s2'),#model4
                    c('sigma1','sigma2','m01','m02','a1','a2','b1','b2','theta1','theta2','alpha1','alpha2','alpha_s1','alpha_s2'),#model5
                    c('sigma1','sigma2','a01','a02','ap1','ap2','lambda1','lambda2','k1','k2','a1','a2','b1','b2','theta1','theta2','alpha1','alpha2','alpha_s1','alpha_s2'),#model6
                    c('sigma1','sigma2','a01','a02','ap1','ap2','lambda1','lambda2','k1','k2','r_win','t_win','alpha1','alpha2','c1','c2','alpha_s1','alpha_s2'),#model7
                    c('sigma1','sigma2','m0_1','m0_2','alpha_c','alpha_m1','alpha_m2','c1','c2','v0_1','v0_2','alpha_s1','alpha_s2'),#model8
                    c('sigma1','sigma2','a01','a02','ap1','ap2','lambda1','lambda2','k1','k2','alpha_c','alpha_m1','alpha_m2','c1','c2','v0_1','v0_2','alpha_s1','alpha_s2'),#model9
                    c('sigma1','sigma2','a01','a02','ap1','ap2','lambda1','lambda2','k1','k2','alpha_c','alpha_m1','alpha_m2','c1','c2','v0_1','v0_2','alpha_s1','alpha_s2','rho_win'),#model10
                    c('NA'),#model11
                    c('NA')#model12
                    )
  
  #range of parameters for each model
  lowerLimit = list(0,#model0
                    c(0,0, 10,10,-100,-100,.1,.1,1,1),#model1
                    c(5,5,20,20,-5,-5,0,0,-1,-1),#model2
                    c(0,0, 10,10,-100,-100,.1,.1,1,1,  -5,-5,0,0,-1,-1),#model3
                    c(0,0, 10,10,-100,-100,.1,.1,1,1,  -5,-5,0,0,-5,-5,-1,-1),#model4
                    c(5,5,20,20,0,0,0,0,0,0,-10,-10,-1,-1),#model5
                    c(0,0, 10,10,-100,-100,.1,.1,1,1,  0,0,0,0,0,0,-10,-10,-1,-1),#model6
                    c(0,0, 10,10,-100,-100,.1,.1,1,1, 1,1,-10,-10,-10,-10, -1,-1),#model7
                    c(0,  0,  0,  0, -.01,-.1,-.1,-20,-20,-5,-5,-1,-1),#model8
                    c(0,0, 10,10,-100,-100,.1,.1,1,1, -.01,-.1,-.1,-20,-20,-5,-5,-1,-1),#model9
                    c(0,0, 10,10,-100,-100,.1,.1,1,1, -.01,-.1,-.1,-20,-20,-5,-5,-1,-1,1),#model10
                    0,#model11
                    0#model12
                    )
  
  upperLimit = list(0,#model0
                    c(100,100,300,300,200,200,20,20,20,20),#model1
                    c(100,100,400,400,15,15,15,5,1,1),#model2
                    c(100,100,300,300,200,200,20,20,20,20,  15,15,15,5,1,1),#model3
                    c(100,100,300,300,200,200,20,20,20,20,  15,15,15,5,5,5,1,1),#model4
                    c(100,100,400,400,20,20,20,20,1,1,10,10,1,1),#model5
                    c(100,100,300,300,200,200,20,20,20,20,  20,20,20,20,1,1,10,10,1,1),#model6
                    c(100,100,300,300,200,200,20,20,20,20,  50,50,10,10,10,10, 1,1),#model7
                    c(100,100,300,300,.01, .1, .1, 20, 20, 5, 5,1,1),#model8
                    c(100,100,300,300,200,200,20,20,20,20,  .01, .1, .1, 20, 20, 5, 5,1,1),#model9
                    c(100,100,300,300,200,200,20,20,20,20,  .01, .1, .1, 20, 20, 5, 5,1,1,50),#model11
                    0,#model11
                    0#model12
                    )
  
}else{#if numOfThreshold==1
  numOfFreeParams = c(2,5,4,0,0,0,0,0,0,0,0,9,12)  
  
  paramNames = list(c('sigma','m0'),#model0
                    c('sigma','a0','ap','lambda','k'),#model1
                    c('sigma','m0','beta','alpha'),#model2
                    c('NA'),#model3
                    c('NA'),#model4
                    c('NA'),#model5
                    c('NA'),#model6
                    c('NA'),#model7
                    c('NA'),#model8
                    c('NA'),#model9
                    c('NA'),#model10
                    c('sigma','a0','ap','lambda','k','r_p','r_n','c','alpha_s'),#model11
                    c('sigma','a0','ap','lambda','k','alpha_c','alpha_m1','alpha_m2','c1','c2','v0_1','v0_2')#model12
  )
  
  #range of parameters for each model
  lowerLimit = list(0,#model0
                    c(0, 10,-100,.1,1),#model1
                    c(5,20,-5,0),#model2
                    0,#model3
                    0,#model4
                    0,#model5
                    0,#model6
                    0,#model7
                    0,#model8
                    0,#model9
                    0,#model10
                    c(0, 10,-100,.1,1,  0,-20,-.1,-1),#model11
                    c(0, 10,-100,.1,1,  -.01,-.1,-.1,-20,-20,-5,-5)#model12
  )
  
  upperLimit = list(0,#model0
                    c(100,300,200,20,20),#model1
                    c(100,400,15,5),#model2
                    0,#model3
                    0,#model4
                    0,#model5
                    0,#model6
                    0,#model7
                    0,#model8
                    0,#model9
                    0,#model10
                    c(100,300,200,20,20, 20,0,.1,1),#model11
                    c(100,300,200,20,20, .01, .1, .1, 20, 20, 5, 5)#model12
                    )
}


















