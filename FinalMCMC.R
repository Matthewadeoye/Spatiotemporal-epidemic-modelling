library(mvtnorm)
library(matrixStats)
library(sf)
library(spdep)
MeningococcalData<- read.csv("Dataset.csv")
head(MeningococcalData)

y<- matrix(MeningococcalData[, 3], nrow = 94, ncol = 156, byrow = T)
e_it<- matrix(MeningococcalData[, 4], nrow = 94, ncol = 156, byrow = T)

France.shp <- st_read("ModifiedFranceShapefile.shp")
France_adjmat <- nb2mat(poly2nb(France.shp), style = "B", zero.policy = TRUE)
init_density<- c(0.6666667, 0.3333333)
time<- ncol(y)
ndept<- nrow(y)
nstate<- 2
R<- -1*France_adjmat
diag(R)<- -rowSums(R, na.rm = T)
qr(R)$rank
source("GeneralLoglikelihood.R")

#Prior density for Trend components (r_t)
randomwalk2<- function(componentR, PrecisionR){
  time<- length(componentR)
  Sumres<- 0
  for(i in 3:time){
    res<- (componentR[i-2] - (2 * componentR[i-1]) + componentR[i])^2
    Sumres<- Sumres + res
  }
  return((time - 2)/2 * log(PrecisionR) - PrecisionR/2 * Sumres)
}

#Prior density for Seasonal components (s_t)
seasonalComp<- function(x, z){
  time<- length(x)
  Sumres<- 0
  for(i in 12:time){
    res<- (sum(x[(i-11):(i-0)]))^2
    Sumres<- Sumres + res
  }
  return((time - 11)/2 * log(z) - z/2 * Sumres)
}

#Intrinsic GMRF density for spatial components (u_i)
logIGMRF1<- function(x, y, z) {
  n = nrow(z)
  sumC = sum(x[1:(n-1)])
  x = c(x[1:(n-1)], -sumC)
  return ((n - 1)/2 * (log(y) - log(2 * pi)) - y/2 * t(x) %*% z %*% x)
}

#Function for transition probability matrix
G<- function(G12, G21){
  m<- matrix(c(1-G12,G12,G21,1-G21),2,2,byrow=T)
  return(m)
}

Model<- 0

num_iteration<- 125000
MC_chain<- matrix(NA, nrow=num_iteration, ncol=413)
MC_chain[1,]<- c(runif(1), runif(1), 8214, 121, 1.4, rep(0, time), rep(0, time), rep(0, ndept), c(0, 0))
acceptedR<- 0
acceptedS<- 0
accepted<- 0

zigmaR<- diag(rep(0.13, time), nrow = time, ncol = time)
zigmaS<- diag(rep(0.13, time), nrow = time, ncol = time)
zigmaU<- diag(rep(0.1, ndept-1), nrow=ndept-1, ncol=ndept-1)
optconstantR<- 2.38^2/(time-2)
optconstantS<- 2.38^2/(time-11)
optconstantU<- 2.38^2/(ndept-1)
lambdaR<- 1
lambdaS<- 1
lambdaU<- 1
lambdakr<- 1
lambdaks<- 1
lambdaku<- 1
sigmakr<- 130000
sigmaks<- 3500
sigmaku<- 2.4

RW2PrecMat<- matrix(0, nrow=156, ncol=156)
RW2PrecMat[1,(1:3)]<- c(1,-2,1)
RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
RW2PrecMat[155,(153:156)]<- c(1,-4,5,-2)
RW2PrecMat[156,(154:156)]<- c(1,-2,1)
for(i in 3:(time-3)){
  RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
}
strr<- RW2PrecMat

A<- 6:18
B<- 19:31
C<- 32:44
D<- 45:57
E<- 58:70
f<- 71:83
g<- 84:96
H<- 97:109
I<- 110:122
J<- 123:135
K<- 136:148
L<- 149:161

Blocks<- list(A,B,C,D,E,f,g,H,I,J,K,L)

SRWPrecMat<- matrix(0, nrow=156, ncol=156)
SRWPrecMat[1,(1:12)]<- rep(1, 12)
SRWPrecMat[2,(1:13)]<- c(1,rep(2, 11),1)
SRWPrecMat[3,(1:14)]<- c(1,2,rep(3, 10),2,1)
SRWPrecMat[4,(1:15)]<- c(1,2,3,rep(4, 9),3,2,1)
SRWPrecMat[5,(1:16)]<- c(1,2,3,4,rep(5, 8),4,3,2,1)
SRWPrecMat[6,(1:17)]<- c(1,2,3,4,5,rep(6, 7),5,4,3,2,1)
SRWPrecMat[7,(1:18)]<- c(1,2,3,4,5,6,rep(7, 6),6,5,4,3,2,1)
SRWPrecMat[8,(1:19)]<- c(1,2,3,4,5,6,7,rep(8, 5),7,6,5,4,3,2,1)
SRWPrecMat[9,(1:20)]<- c(1,2,3,4,5,6,7,8,rep(9, 4),8,7,6,5,4,3,2,1)
SRWPrecMat[10,(1:21)]<- c(1,2,3,4,5,6,7,8,9,rep(10, 3),9,8,7,6,5,4,3,2,1)
SRWPrecMat[11,(1:22)]<- c(1,2,3,4,5,6,7,8,9,10,rep(11, 2),10,9,8,7,6,5,4,3,2,1)
SRWPrecMat[12,(1:23)]<- c(1,2,3,4,5,6,7,8,9,10,11,12,11,10,9,8,7,6,5,4,3,2,1)

for(i in 14:(time-11)){
  SRWPrecMat[i-1, ((i-12):(i+10))]<- c(1,2,3,4,5,6,7,8,9,10,11,12,11,10,9,8,7,6,5,4,3,2,1)
}

SRWPrecMat[145, (134:156)]<-   c(1,2,3,4,5,6,7,8,9,10,11,12,11,10,9,8,7,6,5,4,3,2,1)
SRWPrecMat[146, (135:156)]<-   c(1,2,3,4,5,6,7,8,9,10,rep(11, 2),10,9,8,7,6,5,4,3,2,1)
SRWPrecMat[147, (136:156)]<-   c(1,2,3,4,5,6,7,8,9,rep(10, 3),9,8,7,6,5,4,3,2,1)
SRWPrecMat[148, (137:156)]<-   c(1,2,3,4,5,6,7,8,rep(9, 4),8,7,6,5,4,3,2,1)
SRWPrecMat[149, (138:156)]<-   c(1,2,3,4,5,6,7,rep(8, 5),7,6,5,4,3,2,1)
SRWPrecMat[150, (139:156)]<-   c(1,2,3,4,5,6,rep(7, 6),6,5,4,3,2,1)
SRWPrecMat[151, (140:156)]<-  c(1,2,3,4,5,rep(6, 7),5,4,3,2,1)
SRWPrecMat[152, (141:156)]<- c(1,2,3,4,rep(5, 8),4,3,2,1)
SRWPrecMat[153, (142:156)]<- c(1,2,3,rep(4, 9),3,2,1)
SRWPrecMat[154, (143:156)]<- c(1,2,rep(3, 10),2,1)
SRWPrecMat[155, (144:156)]<- c(1,rep(2, 11),1)
SRWPrecMat[156, (145:156)]<- rep(1, 12)
strs<- SRWPrecMat

SA<- 162:171
SB<- 172:181
SC<- 182:191
SD<- 192:201
SE<- 202:211
Sf<- 212:221
Sg<- 222:231
SH<- 232:241
SI<- 242:251
SJ<- 252:261
SK<- 262:271
SL<- 272:281
SM<- 282:291
SN<- 292:301
SO<- 302:311
SP<- 312:317

SBlocks<- list(SA,SB,SC,SD,SE,Sf,Sg,SH,SI,SJ,SK,SL,SM,SN,SO,SP)

for(i in 2:num_iteration){
  
  if(Model == 0){
    proposedB <- c(0, 0)
  }else if(Model == 1 || Model == 2 || Model == 4 || Model == 5) {
    proposedB <- rnorm(1, mean = MC_chain[i-1, 412], sd = 0.9)
    proposedB <- c(proposedB, 0)
  }else if(Model == 3 || Model == 6){
    proposedB <- rnorm(2, mean = MC_chain[i-1, 412:413], sd = c(0.9, 0.8))
  }
  
  proposedkappaR<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, 6:161]) %*% strr %*% MC_chain[i-1, 6:161])/2)
  
  proposalcurrentkappaR<- dgamma(MC_chain[i-1,3], shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, 6:161]) %*% strr %*% MC_chain[i-1, 6:161])/2, log=TRUE)
  proposalproposedkappaR<- dgamma(proposedkappaR, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, 6:161]) %*% strr %*% MC_chain[i-1, 6:161])/2, log=TRUE)
  
  mh.ratiokr<- exp(proposalproposedkappaR - proposalcurrentkappaR)
  
  print(mh.ratiokr)
  
  if(!is.na(mh.ratiokr) && runif(1) < mh.ratiokr){
    MC_chain[i,3]<- proposedkappaR
    MC_chain[i, 412:413]<- proposedB
    
    accepted<- accepted + 1
  }
  else{
    MC_chain[i,3]<- MC_chain[i-1,3]
    MC_chain[i,412:413]<- MC_chain[i-1,412:413]
    
    accepted<- accepted + 0
  }
  
  proposedkappaS<- rgamma(1, shape = 1 + (time-11)/2, rate = 0.0001 + (t(MC_chain[i-1, 162:317]) %*% strs %*% MC_chain[i-1, 162:317])/2)
  
  proposalcurrentkappas<- dgamma(MC_chain[i-1,4], shape = 1 + (time-11)/2, rate = 0.0001 + (t(MC_chain[i-1, 162:317]) %*% strs %*% MC_chain[i-1, 162:317])/2, log=TRUE)
  proposalproposedkappas<- dgamma(proposedkappaS, shape = 1 + (time-11)/2, rate = 0.0001 + (t(MC_chain[i-1, 162:317]) %*% strs %*% MC_chain[i-1, 162:317])/2, log=TRUE)
  
  mh.ratioks<- exp(proposalproposedkappas - proposalcurrentkappas)

  print(mh.ratioks)
  
  if(!is.na(mh.ratioks) && runif(1) < mh.ratioks){
    MC_chain[i,4]<- proposedkappaS
    
    accepted<- accepted + 1
  }
  else{
    MC_chain[i,4]<- MC_chain[i-1,4]
    
    accepted<- accepted + 0
  }
  
  proposedkappaU<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, 318:411]) %*% R %*% MC_chain[i-1, 318:411])/2)
  
  proposalcurrentkappau<- dgamma(MC_chain[i-1,5], shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, 318:411]) %*% R %*% MC_chain[i-1, 318:411])/2, log=TRUE)
  proposalproposedkappau<- dgamma(proposedkappaU, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, 318:411]) %*% R %*% MC_chain[i-1, 318:411])/2, log=TRUE)
  
  mh.ratioku<- exp(proposalproposedkappau - proposalcurrentkappau)

  print(mh.ratioku)
  
  if(!is.na(mh.ratioku) && runif(1) < mh.ratioku){
    MC_chain[i,5]<- proposedkappaU
    
    accepted<- accepted + 1
  }
  else{
    MC_chain[i,5]<- MC_chain[i-1,5]
    
    accepted<- accepted + 0
  }
  
  proposedG12<- abs(rnorm(1,mean=MC_chain[i-1,1], sd=0.9))
  if(proposedG12>1) proposedG12=2-proposedG12
  
  priorcurrent12<- dbeta(MC_chain[i-1,1], shape1 = 1, shape2 = 1, log=TRUE)
  priorproposed12<- dbeta(proposedG12, shape1 = 1, shape2 = 1, log=TRUE) 
  
  likelihoodcurrent<- GeneralLoglikelihood(y, MC_chain[i-1,6:161], MC_chain[i-1,162:317], MC_chain[i-1, 318:411], G(MC_chain[i-1,1],MC_chain[i-1,2]),init_density,e_it, MC_chain[i-1, 412:413], Model, France_adjmat)
  likelihoodproposed<- GeneralLoglikelihood(y, MC_chain[i-1, 6:161], MC_chain[i-1, 162:317], MC_chain[i-1, 318:411], G(proposedG12,MC_chain[i-1,2]),init_density,e_it, MC_chain[i-1, 412:413], Model,France_adjmat)
  
  mh.ratio<- exp(likelihoodproposed + priorproposed12
                 - likelihoodcurrent - priorcurrent12)
  
  print(mh.ratio)
  
  if(!is.na(mh.ratio) && runif(1) < mh.ratio){
    MC_chain[i, 1]<- proposedG12
    
    accepted<- accepted + 1
  }
  else{
    MC_chain[i, 1]<- MC_chain[i-1,1]
    
    accepted<- accepted + 0
  }
  
  proposedG21<- abs(rnorm(1,mean=MC_chain[i-1,2], sd=0.9))
  if(proposedG21>1) proposedG21=2-proposedG21
  
  priorcurrent21<- dbeta(MC_chain[i-1,2], shape1 = 1, shape2 = 1, log=TRUE)
  priorproposed21<- dbeta(proposedG21, shape1 = 1, shape2 = 1, log=TRUE)
  
  likelihoodcurrent<- GeneralLoglikelihood(y, MC_chain[i-1,6:161], MC_chain[i-1,162:317], MC_chain[i-1, 318:411], G(MC_chain[i,1],MC_chain[i-1,2]),init_density,e_it, MC_chain[i-1, 412:413], Model, France_adjmat)
  likelihoodproposed<- GeneralLoglikelihood(y, MC_chain[i-1, 6:161], MC_chain[i-1, 162:317], MC_chain[i-1, 318:411], G(MC_chain[i,1],proposedG21),init_density,e_it, MC_chain[i-1, 412:413], Model,France_adjmat)
  
  mh.ratio<- exp(likelihoodproposed + priorproposed21
                 - likelihoodcurrent - priorcurrent21)
  
  print(mh.ratio)
  
  if(!is.na(mh.ratio) && runif(1) < mh.ratio){
    MC_chain[i,2]<- proposedG21
    
    accepted<- accepted + 1
  }
  else{
    MC_chain[i,2]<- MC_chain[i-1,2]
    
    accepted<- accepted + 0
  }
  
  proposedspatcomps<- rmvnorm(1, mean=MC_chain[i-1, 318:410], sigma = zigmaU) 
  proposedspatcomps<- c(proposedspatcomps, -sum(proposedspatcomps))
  
  priorcurrentUcomps<- logIGMRF1(MC_chain[i-1, 318:411], MC_chain[i, 5], R)
  priorproposedUcomps<- logIGMRF1(proposedspatcomps, MC_chain[i, 5], R)
  
  proposalproposedcompsU<- sum(dmvnorm(proposedspatcomps[-94], mean = MC_chain[i-1, 318:410], sigma = zigmaU, log = TRUE))
  proposalcurrentcompsU<- sum(dmvnorm(MC_chain[i-1, 318:410], mean = proposedspatcomps[-94], sigma = zigmaU, log = TRUE))
  
  likelihoodcurrent<- GeneralLoglikelihood(y,MC_chain[i-1,6:161],MC_chain[i-1,162:317],MC_chain[i-1,318:411],G(MC_chain[i,1],MC_chain[i,2]),init_density,e_it, MC_chain[i-1, 412:413], Model,France_adjmat)
  likelihoodproposed<- GeneralLoglikelihood(y,MC_chain[i-1,6:161],MC_chain[i-1,162:317],proposedspatcomps,G(MC_chain[i,1],MC_chain[i,2]),init_density,e_it, MC_chain[i-1, 412:413], Model,France_adjmat)
  
  mh.ratioU<- exp(likelihoodproposed + priorproposedUcomps + proposalcurrentcompsU
                  - likelihoodcurrent - priorcurrentUcomps - proposalproposedcompsU)
  
  print(mh.ratioU)
  
  if(!is.na(mh.ratioU) && runif(1) < mh.ratioU){
    MC_chain[i,318:411]<- proposedspatcomps
    
    accepted<- accepted + 1
  }
  else{
    MC_chain[i,318:411]<- MC_chain[i-1,318:411]
    
    accepted<- accepted + 0
  }
  
  RW2PrecMat<- MC_chain[i, 3] * strr
  RconditionalcovA<- solve(RW2PrecMat[1:13, 1:13]) 
  RconditionalcovB<- solve(RW2PrecMat[14:26, 14:26]) 
  RconditionalcovC<- solve(RW2PrecMat[27:39, 27:39])
  RconditionalcovD<- solve(RW2PrecMat[40:52, 40:52])  
  RconditionalcovE<- solve(RW2PrecMat[53:65, 53:65]) 
  RconditionalcovF<- solve(RW2PrecMat[66:78, 66:78])
  RconditionalcovG<- solve(RW2PrecMat[79:91, 79:91]) 
  RconditionalcovH<- solve(RW2PrecMat[92:104, 92:104]) 
  RconditionalcovI<- solve(RW2PrecMat[105:117, 105:117]) 
  RconditionalcovJ<- solve(RW2PrecMat[118:130, 118:130]) 
  RconditionalcovK<- solve(RW2PrecMat[131:143, 131:143]) 
  RconditionalcovL<- solve(RW2PrecMat[144:156, 144:156]) 
  
  covBlocks<- list(RconditionalcovA, RconditionalcovB, RconditionalcovC, RconditionalcovD, RconditionalcovE, RconditionalcovF, RconditionalcovG,
                   RconditionalcovH, RconditionalcovI, RconditionalcovJ, RconditionalcovK, RconditionalcovL)
  
  
  SRWPrecMat<- MC_chain[i, 4] * strs
  SconditionalcovA<- solve(SRWPrecMat[1:10, 1:10]) 
  SconditionalcovB<- solve(SRWPrecMat[11:20, 11:20]) 
  SconditionalcovC<- solve(SRWPrecMat[21:30, 21:30])
  SconditionalcovD<- solve(SRWPrecMat[31:40, 31:40])  
  SconditionalcovE<- solve(SRWPrecMat[41:50, 41:50]) 
  SconditionalcovF<- solve(SRWPrecMat[51:60, 51:60])
  SconditionalcovG<- solve(SRWPrecMat[61:70, 61:70]) 
  SconditionalcovH<- solve(SRWPrecMat[71:80, 71:80]) 
  SconditionalcovI<- solve(SRWPrecMat[81:90, 81:90]) 
  SconditionalcovJ<- solve(SRWPrecMat[91:100, 91:100]) 
  SconditionalcovK<- solve(SRWPrecMat[101:110, 101:110]) 
  SconditionalcovL<- solve(SRWPrecMat[111:120, 111:120]) 
  SconditionalcovM<- solve(SRWPrecMat[121:130, 121:130]) 
  SconditionalcovN<- solve(SRWPrecMat[131:140, 131:140]) 
  SconditionalcovO<- solve(SRWPrecMat[141:150, 141:150]) 
  SconditionalcovP<- solve(SRWPrecMat[151:156, 151:156]) 
  
  ScovBlocks<- list(SconditionalcovA, SconditionalcovB, SconditionalcovC, SconditionalcovD, SconditionalcovE, SconditionalcovF, SconditionalcovG,
                    SconditionalcovH, SconditionalcovI, SconditionalcovJ, SconditionalcovK, SconditionalcovL, SconditionalcovM, SconditionalcovN, SconditionalcovO, SconditionalcovP)
  
  for(j in 1:12){
    if(j==1){
      Rconditionalmean<- -covBlocks[[j]] %*% RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[-j]))-5] %*% MC_chain[i-1, unlist(Blocks[-j])]
      proposedRcomps<- rmvnorm(1, mean = Rconditionalmean, sigma = covBlocks[[j]]) 
      proposedRcomps<- c(proposedRcomps, MC_chain[i-1, unlist(Blocks[-j])])
      
    }
    else if(j!=1 && j!=12){ 
      Rconditionalmean<- -covBlocks[[j]] %*% (RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[1:(j-1)]))-5] %*% MC_chain[i, unlist(Blocks[1:(j-1)])] + RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[(j+1):12]))-5] %*% MC_chain[i, unlist(Blocks[(j+1):12])]) 
      proposedRcomps<- rmvnorm(1, mean = Rconditionalmean, sigma = covBlocks[[j]]) 
      proposedRcomps<- c(MC_chain[i, unlist(Blocks[1:(j-1)])], proposedRcomps, MC_chain[i, unlist(Blocks[(j+1):12])])
      
    }
    else if(j==12){
      Rconditionalmean<- -covBlocks[[j]] %*% RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[-j]))-5] %*% MC_chain[i, unlist(Blocks[-j])]
      proposedRcomps<- rmvnorm(1, mean = Rconditionalmean, sigma = covBlocks[[j]]) 
      proposedRcomps<- c(MC_chain[i, unlist(Blocks[-j])], proposedRcomps)
      
    }
    
    likelihoodcurrent<- GeneralLoglikelihood(y,MC_chain[i-1, 6:161],MC_chain[i-1,162:317],MC_chain[i,318:411],G(MC_chain[i,1],MC_chain[i,2]),init_density,e_it, MC_chain[i-1, 412:413], Model,France_adjmat)
    likelihoodproposed<- GeneralLoglikelihood(y,proposedRcomps,MC_chain[i-1,162:317],MC_chain[i,318:411],G(MC_chain[i,1],MC_chain[i,2]),init_density,e_it, MC_chain[i-1, 412:413], Model,France_adjmat)
    
    mh.ratioR<- exp(likelihoodproposed - likelihoodcurrent)
    
    print(mh.ratioR)
    if(!is.na(mh.ratioR) && runif(1) < mh.ratioR){
      MC_chain[i,6:161]<- proposedRcomps

    }
    else{
      if(j==1){
        MC_chain[i,6:161]<- MC_chain[i-1,6:161]
        
      }
      else if(j!=1){
        MC_chain[i,6:161]<- MC_chain[i,6:161]
      }      
    }
  }
  
  for(j in 1:16){
    if(j==1){
      Sconditionalmean<- -ScovBlocks[[j]] %*% SRWPrecMat[(SBlocks[[j]])-161, (unlist(SBlocks[-j]))-161] %*% MC_chain[i-1, unlist(SBlocks[-j])]
      proposedScomps<- rmvnorm(1, mean = Sconditionalmean, sigma = ScovBlocks[[j]]) 
      proposedScomps<- c(proposedScomps, MC_chain[i-1, unlist(SBlocks[-j])])
    }
    else if(j!=1 && j!=16){ 
      Sconditionalmean<- -ScovBlocks[[j]] %*% (SRWPrecMat[(SBlocks[[j]])-161, (unlist(SBlocks[1:(j-1)]))-161] %*% MC_chain[i, unlist(SBlocks[1:(j-1)])] + SRWPrecMat[(SBlocks[[j]])-161, (unlist(SBlocks[(j+1):16]))-161] %*% MC_chain[i, unlist(SBlocks[(j+1):16])]) 
      proposedScomps<- rmvnorm(1, mean = Sconditionalmean, sigma = ScovBlocks[[j]]) 
      proposedScomps<- c(MC_chain[i, unlist(SBlocks[1:(j-1)])], proposedScomps, MC_chain[i, unlist(SBlocks[(j+1):16])])
    }
    else if(j==16){
      Sconditionalmean<- -ScovBlocks[[j]] %*% SRWPrecMat[(SBlocks[[j]])-161, (unlist(SBlocks[-j]))-161] %*% MC_chain[i, unlist(SBlocks[-j])]
      proposedScomps<- rmvnorm(1, mean = Sconditionalmean, sigma = ScovBlocks[[j]]) 
      proposedScomps<- c(MC_chain[i, unlist(SBlocks[-j])], proposedScomps)
    }
    
    likelihoodcurrent<- GeneralLoglikelihood(y,MC_chain[i, 6:161],MC_chain[i-1,162:317],MC_chain[i,318:411],G(MC_chain[i,1],MC_chain[i,2]),init_density,e_it, MC_chain[i-1, 412:413], Model,France_adjmat)
    likelihoodproposed<- GeneralLoglikelihood(y,MC_chain[i, 6:161],proposedScomps,MC_chain[i,318:411],G(MC_chain[i,1],MC_chain[i,2]),init_density,e_it, MC_chain[i-1, 412:413], Model,France_adjmat)
    
    mh.ratioS<- exp(likelihoodproposed - likelihoodcurrent)
    
    print(mh.ratioS)
    if(!is.na(mh.ratioS) && runif(1) < mh.ratioS){
      MC_chain[i,162:317]<- proposedScomps
    }
    else{
      if(j==1){
        MC_chain[i,162:317]<- MC_chain[i-1,162:317]
        
      }
      else if(j!=1){
        MC_chain[i,162:317]<- MC_chain[i,162:317]
      }      
    }
  }
  
  #Adapting zigmaR
  if(i==10000){
    epsilonR<- 0.7
    XnR<- MC_chain[1:i, 6:161]
    XnbarR <- colMeans(XnR) 
    zigmaR <- cov(XnR) + epsilonR * diag(rep(1, time))
    zigmaR<- optconstantR * zigmaR
  } else if (i > 10000){ 
    
    ### Using random walk after 10000 conditional prior proposals
    
    proposedRcomps<- rmvnorm(1, mean = MC_chain[i, 6:161], sigma = zigmaR)
    
    priorcurrentRcomps<- randomwalk2(MC_chain[i, 6:161], MC_chain[i, 3])
    priorproposedRcomps<- randomwalk2(proposedRcomps, MC_chain[i, 3])
    
    proposalproposedRcomps<- sum(dmvnorm(proposedRcomps, mean = MC_chain[i, 6:161], sigma = zigmaR, log = TRUE))
    proposalcurrentRcomps<- sum(dmvnorm(MC_chain[i, 6:161], mean = proposedRcomps, sigma = zigmaR, log = TRUE))
    
    likelihoodcurrent<- GeneralLoglikelihood(y,MC_chain[i,6:161],MC_chain[i, 162:317],MC_chain[i,318:411],G(MC_chain[i,1],MC_chain[i,2]),init_density,e_it, MC_chain[i-1, 412:413], Model,France_adjmat)
    likelihoodproposed<- GeneralLoglikelihood(y,proposedRcomps,MC_chain[i, 162:317],MC_chain[i,318:411],G(MC_chain[i,1],MC_chain[i,2]),init_density,e_it, MC_chain[i-1, 412:413], Model,France_adjmat)
    
    mh.ratioR<- exp(likelihoodproposed + priorproposedRcomps + proposalcurrentRcomps
                    - likelihoodcurrent - priorcurrentRcomps - proposalproposedRcomps)
    
    print(mh.ratioR)
    
    if(!is.na(mh.ratioR) && runif(1) < mh.ratioR){
      MC_chain[i,6:161]<- proposedRcomps
      
      acceptedR<- acceptedR + 1
    }
    else{
      MC_chain[i,6:161]<- MC_chain[i,6:161]
      
      acceptedR<- acceptedR + 0
    }
    
    proposedScomps<- rmvnorm(1, mean = MC_chain[i, 162:317], sigma = zigmaS)
    
    priorcurrentScomps<- seasonalComp(MC_chain[i, 162:317], MC_chain[i, 4])
    priorproposedScomps<- seasonalComp(proposedScomps, MC_chain[i, 4])
    
    proposalproposedScomps<- sum(dmvnorm(proposedScomps, mean = MC_chain[i, 162:317], sigma = zigmaS, log = TRUE))
    proposalcurrentScomps<- sum(dmvnorm(MC_chain[i, 162:317], mean = proposedScomps, sigma = zigmaS, log = TRUE))
    
    likelihoodcurrent<- GeneralLoglikelihood(y,MC_chain[i,6:161],MC_chain[i,162:317],MC_chain[i,318:411],G(MC_chain[i,1], MC_chain[i,2]),init_density,e_it, MC_chain[i-1, 412:413], Model,France_adjmat)
    likelihoodproposed<- GeneralLoglikelihood(y,MC_chain[i,6:161],proposedScomps,MC_chain[i,318:411],G(MC_chain[i,1], MC_chain[i,2]),init_density,e_it, MC_chain[i-1, 412:413], Model,France_adjmat)
    
    mh.ratioS<- exp(likelihoodproposed + priorproposedScomps + proposalcurrentScomps
                    - likelihoodcurrent - priorcurrentScomps - proposalproposedScomps)
    
    print(mh.ratioS)
    
    if(!is.na(mh.ratioS) && runif(1) < mh.ratioS){
      MC_chain[i,162:317]<- proposedScomps
      
      acceptedS<- acceptedS + 1
    }
    else{
      MC_chain[i,162:317]<- MC_chain[i,162:317]
      
      acceptedS<- acceptedS + 0
    }
    
    XnbarPrevR <- XnbarR
    XnbarR <- (i*XnbarR + MC_chain[i, 6:161])/(i+1)
    zigmaR <- ((i-1)*zigmaR + tcrossprod(MC_chain[i, 6:161]) + i*tcrossprod(XnbarPrevR) - (i+1)*tcrossprod(XnbarR) + epsilonR*diag(rep(1,time)))/i
    #Robbins Munro tuning
    lambdaR<- lambdaR * exp((2/max(1, i-10000)) * (min(mh.ratioR, 1) - 0.234))
    zigmaR<- lambdaR* optconstantR * zigmaR
    #print(zigmaR)
  }   
  
  #Adapting zigmaS
  if(i==10000){
    epsilonS<- 0.7
    XnS<- MC_chain[1:i, 162:317]
    XnbarS <- colMeans(XnS) 
    zigmaS <- cov(XnS) + epsilonS*diag(rep(1, time))
    zigmaS<- optconstantS * zigmaS
  } else if (i > 10000){ 
    XnbarPrevS <- XnbarS
    XnbarS <- (i*XnbarS + MC_chain[i, 162:317])/(i+1)
    zigmaS <- ((i-1)*zigmaS + tcrossprod(MC_chain[i, 162:317]) + i*tcrossprod(XnbarPrevS) - (i+1)*tcrossprod(XnbarS) + epsilonS*diag(rep(1,time)))/i
    #Robbins Munro tuning
    lambdaS<- lambdaS * exp((2/max(1, i-10000)) * (min(mh.ratioS, 1) - 0.234))
    zigmaS<- lambdaS* optconstantS * zigmaS
    #print(zigmaS)
  }   
  
  #Adapting zigmaU
  if(i==10000){
    epsilonU<- 0.7
    XnU<- MC_chain[1:i, 318:410]
    XnbarU <- colMeans(XnU) 
    zigmaU <- cov(XnU) + epsilonU*diag(rep(1, ndept-1))
    zigmaU<- optconstantU * zigmaU
  } else if (i > 10000){ 
    XnbarPrevU <- XnbarU
    XnbarU <- (i*XnbarU + MC_chain[i, 318:410])/(i+1)
    zigmaU <- ((i-1)*zigmaU + tcrossprod(MC_chain[i, 318:410]) + i*tcrossprod(XnbarPrevU) - (i+1)*tcrossprod(XnbarU) + epsilonU*diag(rep(1,ndept-1)))/i
    #Robbins Munro tuning
    lambdaU<- lambdaU * exp((2/max(1, i-10000)) * (min(mh.ratioU, 1) - 0.234))
    zigmaU<- lambdaU* optconstantU * zigmaU
  } 
}

#Results
RsMCMC<- numeric(time)
for(i in 1:time){
  RsMCMC[i] = mean(MC_chain[,5+i])
}
RsMCMC

SsMCMC<- numeric(time)
for(i in 1:time){
  SsMCMC[i] = mean(MC_chain[,161+i])
}
SsMCMC

UsMCMC<- numeric(ndept)
for(i in 1:ndept){
  UsMCMC[i] = mean(MC_chain[,317+i])
}
UsMCMC

acceptedR/num_iteration
acceptedS/num_iteration
accepted/num_iteration
