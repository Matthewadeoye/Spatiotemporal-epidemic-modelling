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
time<- ncol(y)
ndept<- nrow(y)
nstate<- 2
R<- -1*France_adjmat
diag(R)<- -rowSums(R, na.rm = T)
qr(R)$rank
source("GeneralLoglikelihood.R")
source("DesignMatrix.R")
source("OutbreakProbability.R")

#Crude estimates
ydot<- colSums(y)
edot<- colSums(e_it)
logydot<- log(ydot)
logedot<- log(edot)
lambdadot<- logydot-logedot
x<- 1:time
loess_fit <- loess(lambdadot ~ x, span = 0.5)
smoothed <- predict(loess_fit)
crudeS<- lambdadot - smoothed
crudeR<- smoothed
crudeU<- log(rowSums(y/e_it)/sum(exp(crudeR+crudeS)))-mean(log(rowSums(y/e_it)/sum(exp(crudeR+crudeS))))

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

Model<- 1
z_it<- DesignMatrixModel1(y, France_adjmat)[[1]]
z_it2<- DesignMatrixModel1(y, France_adjmat)[[2]]

num_iteration<- 25000
MC_chain<- matrix(NA, nrow=num_iteration, ncol=5+time+time+ndept+2)
MC_chain[1,]<- c(runif(1), runif(1), 8214, 121, 1.4, crudeR, crudeS, crudeU, c(0, 0))
acceptedR<- 0
acceptedS<- 0
accepted<- 0
#likelihoodcurrent<- GeneralLoglikelihood(y, MC_chain[1,5+(1:time)], MC_chain[1, 5+time+(1:time)], MC_chain[1, 318:411], G(MC_chain[1, 1],MC_chain[1, 2]), e_it, MC_chain[1, 5+time+time+ndept+(1:2)], Model, z_it, z_it2)

zigmaR<- diag(rep(0.13, time), nrow = time, ncol = time)
zigmaS<- diag(rep(0.13, time), nrow = time, ncol = time)
zigmaU<- diag(rep(0.1, ndept-1), nrow=ndept-1, ncol=ndept-1)
optconstantR<- 2.38^2/(time-2)
optconstantS<- 2.38^2/(time-11)
optconstantU<- 2.38^2/(ndept-1)
lambdaR<- 1
lambdaS<- 1
lambdaU<- 1

RW2PrecMat<- matrix(0, nrow=time, ncol=time)
RW2PrecMat[1,(1:3)]<- c(1,-2,1)
RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
for(i in 3:(time-3)){
  RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
}
strr<- RW2PrecMat

A<- 5+(1:13)
B<- 5+(14:26)
C<- 5+(27:39)
D<- 5+(40:52)
E<- 5+(53:65)
f<- 5+(66:78)
g<- 5+(79:91)
H<- 5+(92:104)
I<- 5+(105:117)
J<- 5+(118:130)
K<- 5+(131:143)
L<- 5+(144:156)
Blocks<- list(A,B,C,D,E,f,g,H,I,J,K,L)

SRWPrecMat<- matrix(0, nrow=time, ncol=time)
SRWPrecMat[1,(1:12)]<- rep(1, 12)
for(i in 2:12){
  SRWPrecMat[i, 1:(i+11)]<- c(1:(i-1), rep(i, (12+1-i)), (i-1):1)
}
for(i in 14:(time-11)){
  SRWPrecMat[i-1, ((i-12):(i+10))]<- c(1:12,11:1)
}
for(i in (time-1):(time-11)){
  SRWPrecMat[i, (i-11):time]<- c(1:(time-i), rep((time-i)+1, i-(time-12)), (time-i):1)
}
SRWPrecMat[time, ((time-11):time)]<- rep(1, 12)
strs<- SRWPrecMat

size<- 10
SA<- 5+time+(1:size)
SB<- 5+time+((size+1):(2*size))
SC<- 5+time+((2*size+1):(3*size))
SD<- 5+time+((3*size+1):(4*size))
SE<- 5+time+((4*size+1):(5*size))
Sf<- 5+time+((5*size+1):(6*size))
Sg<- 5+time+((6*size+1):(7*size))
SH<- 5+time+((7*size+1):(8*size))
SI<- 5+time+((8*size+1):(9*size))
SJ<- 5+time+((9*size+1):(10*size))
SK<- 5+time+((10*size+1):(11*size))
SL<- 5+time+((11*size+1):(12*size))
SM<- 5+time+((12*size+1):(13*size))
SN<- 5+time+((13*size+1):(14*size))
SO<- 5+time+((14*size+1):(15*size))
SP<- 5+time+((15*size+1):time)
SBlocks<- list(SA,SB,SC,SD,SE,Sf,Sg,SH,SI,SJ,SK,SL,SM,SN,SO,SP)

for(i in 2:num_iteration){
  print(i)
  
  proposedkappaR<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, 5+(1:time)]) %*% strr %*% MC_chain[i-1, 5+(1:time)])/2)
  MC_chain[i,3]<- proposedkappaR
  print(paste("GibbskappaR = ", proposedkappaR))
  
  proposedkappaS<- rgamma(1, shape = 1 + (time-11)/2, rate = 0.0001 + (t(MC_chain[i-1, 5+time+(1:time)]) %*% strs %*% MC_chain[i-1, 5+time+(1:time)])/2)
  MC_chain[i,4]<- proposedkappaS
  print(paste("GibbskappaS = ", proposedkappaS))
  
  proposedkappaU<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, 5+time+time+(1:ndept)]) %*% R %*% MC_chain[i-1, 5+time+time+(1:ndept)])/2)
  MC_chain[i,5]<- proposedkappaU
  print(paste("GibbskappaU = ", proposedkappaU))
  
  proposedspatcomps<- rmvnorm(1, mean=MC_chain[i-1, 5+time+time+(1:(ndept-1))], sigma = zigmaU) 
  proposedspatcomps<- c(proposedspatcomps, -sum(proposedspatcomps))
  
  priorcurrentUcomps<- logIGMRF1(MC_chain[i-1, 5+time+time+(1:ndept)], MC_chain[i, 5], R)
  priorproposedUcomps<- logIGMRF1(proposedspatcomps, MC_chain[i, 5], R)
  
  proposalproposedcompsU<- sum(dmvnorm(proposedspatcomps[-94], mean = MC_chain[i-1, 5+time+time+(1:(ndept-1))], sigma = zigmaU, log = TRUE))
  proposalcurrentcompsU<- sum(dmvnorm(MC_chain[i-1, 5+time+time+(1:(ndept-1))], mean = proposedspatcomps[-94], sigma = zigmaU, log = TRUE))
  
  likelihoodcurrent<- GeneralLoglikelihood(y,MC_chain[i-1, 5+(1:time)],MC_chain[i-1, 5+time+(1:time)], MC_chain[i-1, 5+time+time+(1:ndept)], G(MC_chain[i-1, 1],MC_chain[i-1,2]), e_it, MC_chain[i-1, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)
  likelihoodproposed<- GeneralLoglikelihood(y,MC_chain[i-1, 5+(1:time)],MC_chain[i-1, 5+time+(1:time)], proposedspatcomps, G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)
  
  mh.ratioU<- exp(likelihoodproposed + priorproposedUcomps + proposalcurrentcompsU
                  - likelihoodcurrent - priorcurrentUcomps - proposalproposedcompsU)
  
  print(paste("mh.ratioU spatcomps =", mh.ratioU))
  
  if(!is.na(mh.ratioU) && runif(1) < mh.ratioU){
    MC_chain[i,5+time+time+(1:ndept)]<- proposedspatcomps
    
    accepted<- accepted + 1
  }
  else{
    MC_chain[i,5+time+time+(1:ndept)]<- MC_chain[i-1,5+time+time+(1:ndept)]
    
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
    
    likelihoodcurrent<- GeneralLoglikelihood(y, MC_chain[i-1, 5+(1:time)], MC_chain[i-1, 5+time+(1:time)], MC_chain[i, 5+time+time+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+time+ndept+(1:2)], Model, z_it, z_it2)
    likelihoodproposed<- GeneralLoglikelihood(y, proposedRcomps, MC_chain[i-1, 5+time+(1:time)], MC_chain[i, 5+time+time+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+time+ndept+(1:2)], Model, z_it, z_it2)
    
    mh.ratioR<- exp(likelihoodproposed - likelihoodcurrent)
    
    print(paste("mh.ratioR condpriorprop =", mh.ratioR))
    if(!is.na(mh.ratioR) && runif(1) < mh.ratioR){
      MC_chain[i, 5+(1:time)]<- proposedRcomps
    }
    else{
      if(j==1){
        MC_chain[i, 5+(1:time)]<- MC_chain[i-1, 5+(1:time)]        
      }
      else if(j!=1){
        MC_chain[i, 5+(1:time)]<- MC_chain[i, 5+(1:time)] 
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
    
    likelihoodcurrent<- GeneralLoglikelihood(y,MC_chain[i, 5+(1:time)],MC_chain[i-1,5+time+(1:time)],MC_chain[i, 5+time+time+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)
    likelihoodproposed<- GeneralLoglikelihood(y,MC_chain[i, 5+(1:time)],proposedScomps,MC_chain[i, 5+time+time+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)
    
    mh.ratioS<- exp(likelihoodproposed - likelihoodcurrent)
    
    print(paste("mh.ratioS condpriorprop = ", mh.ratioS))
    if(!is.na(mh.ratioS) && runif(1) < mh.ratioS){
      MC_chain[i, 5+time+(1:time)]<- proposedScomps
    }
    else{
      if(j==1){
        MC_chain[i, 5+time+(1:time)]<- MC_chain[i-1, 5+time+(1:time)] 
      }
      else if(j!=1){
        MC_chain[i, 5+time+(1:time)]<- MC_chain[i, 5+time+(1:time)]       
      }      
    }
  }
  
  if(Model == 0){
    proposedB <- c(0, 0)
  }else if(Model == 1 || Model == 2 || Model == 4 || Model == 5) {
    proposedB <- abs(rnorm(1, mean = MC_chain[i-1, 412], sd = 0.09))
    proposedB <- c(proposedB, 0)
  }else if(Model == 3 || Model == 6){
    proposedB <- abs(rnorm(2, mean = MC_chain[i-1, 5+time+time+ndept+(1:2)], sd = c(0.09, 0.08)))
  }
  
  priorcurrentB<- sum(dexp(MC_chain[i-1, 5+time+time+ndept+(1:2)], rate = 1, log=TRUE))
  priorproposedB<- sum(dexp(proposedB, rate = 1, log=TRUE)) 
  
  likelihoodcurrent<- GeneralLoglikelihood(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:time)], MC_chain[i, 5+time+time+(1:ndept)], G(MC_chain[i-1, 1],MC_chain[i-1,2]), e_it, MC_chain[i-1, 5+time+time+ndept+(1:2)], Model, z_it, z_it2)
  likelihoodproposed<- GeneralLoglikelihood(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:time)], MC_chain[i, 5+time+time+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1,2]), e_it, proposedB, Model,z_it, z_it2)
  
  mh.ratio<- exp(likelihoodproposed + priorproposedB
                 - likelihoodcurrent - priorcurrentB)
  
  print(paste("mh.ratioB = ", mh.ratio))
  
  if(!is.na(mh.ratio) && runif(1) < mh.ratio){
    MC_chain[i, 5+time+time+ndept+(1:2)]<- proposedB
  }
  else{
    MC_chain[i, 5+time+time+ndept+(1:2)]<- MC_chain[i-1, 5+time+time+ndept+(1:2)]
  }
  
  MC<- Decoding(y = y, r = MC_chain[i, 5+(1:time)], s = MC_chain[i, 5+time+(1:time)], u = MC_chain[i, 5+time+time+(1:ndept)], Gamma = G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it = e_it, B = MC_chain[i, 5+time+time+ndept+(1:2)], Model = Model, z_it = z_it, z_it2 = z_it2)
  ifelse(MC > 0.5, 1, 0)
  MC11<- matrix(0, nrow = ndept, ncol = time)
  MC12<- matrix(0, nrow = ndept, ncol = time)
  MC21<- matrix(0, nrow = ndept, ncol = time)
  MC22<- matrix(0, nrow = ndept, ncol = time)
  
  for(l in 1:ndept){
    for(t in 1:(time-1)){
      if(MC[l, t]== 0 && MC[l, t+1]== 0){
        MC11[l, t]<- 1
      }else{
        MC11[l, t]<- 0
      }
    }
  }
  
  for(l in 1:ndept){
    for(t in 1:(time-1)){
      if(MC[l, t]== 0 && MC[l, t+1]== 1){
        MC12[l, t]<- 1
      }else{
        MC12[l, t]<- 0
      }
    }
  }  
  
  for(l in 1:ndept){
    for(t in 1:(time-1)){
      if(MC[l, t]== 1 && MC[l, t+1]== 0){
        MC21[l, t]<- 1
      }else{
        MC21[l, t]<- 0
      }
    }
  }  
  
  for(l in 1:ndept){
    for(t in 1:(time-1)){
      if(MC[l, t]== 1 && MC[l, t+1]== 1){
        MC22[l, t]<- 1
      }else{
        MC22[l, t]<- 0
      }
    }
  }  
  
  proposedG12<- rbeta(1, shape1 = 1 + sum(MC12), shape2 = 1 + sum(MC11))
  MC_chain[i, 1]<- proposedG12
  print(paste("GibbsG12 = ", proposedG12))
  
  proposedG21<- rbeta(1, shape1 = 1 + sum(MC21), shape2 = 1 + sum(MC22))
  MC_chain[i, 2]<- proposedG21
  print(paste("GibbsG21 = ", proposedG21))
  
  #Adapting zigmaR
  if(i==10){
    epsilonR<- 0.000007
    XnR<- MC_chain[1:i, 5+(1:time)]
    XnbarR <- colMeans(XnR) 
    zigmaR <- cov(XnR) + epsilonR * diag(rep(1, time))
    zigmaR<- optconstantR * zigmaR
  }else if(i > 10){ 
    
    ### Start random walk after 10 conditional prior proposals
    
    proposedRcomps<- rmvnorm(1, mean = MC_chain[i, 5+(1:time)], sigma = zigmaR)
    
    priorcurrentRcomps<- randomwalk2(MC_chain[i, 5+(1:time)], MC_chain[i, 3])
    priorproposedRcomps<- randomwalk2(proposedRcomps, MC_chain[i, 3])
    
    proposalproposedRcomps<- sum(dmvnorm(proposedRcomps, mean = MC_chain[i, 5+(1:time)], sigma = zigmaR, log = TRUE))
    proposalcurrentRcomps<- sum(dmvnorm(MC_chain[i, 5+(1:time)], mean = proposedRcomps, sigma = zigmaR, log = TRUE))
    
    likelihoodcurrent<- GeneralLoglikelihood(y,MC_chain[i, 5+(1:time)],MC_chain[i, 5+time+(1:time)],MC_chain[i, 5+time+time+(1:ndept)],G(MC_chain[i,1],MC_chain[i,2]), e_it, MC_chain[i, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)
    likelihoodproposed<- GeneralLoglikelihood(y,proposedRcomps,MC_chain[i, 5+time+(1:time)],MC_chain[i, 5+time+time+(1:ndept)],G(MC_chain[i,1],MC_chain[i,2]), e_it, MC_chain[i, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)
    
    mh.ratioR<- exp(likelihoodproposed + priorproposedRcomps + proposalcurrentRcomps
                    - likelihoodcurrent - priorcurrentRcomps - proposalproposedRcomps)
    
    print(paste("mh.ratioR = ", mh.ratioR))
    
    if(!is.na(mh.ratioR) && runif(1) < mh.ratioR){
      MC_chain[i, 5+(1:time)]<- proposedRcomps
      
      acceptedR<- acceptedR + 1
    }
    else{
      MC_chain[i, 5+(1:time)]<- MC_chain[i, 5+(1:time)]
      
      acceptedR<- acceptedR + 0
    }
    
    proposedScomps<- rmvnorm(1, mean = MC_chain[i, 5+time+(1:time)], sigma = zigmaS)
    
    priorcurrentScomps<- seasonalComp(MC_chain[i, 5+time+(1:time)], MC_chain[i, 4])
    priorproposedScomps<- seasonalComp(proposedScomps, MC_chain[i, 4])
    
    proposalproposedScomps<- sum(dmvnorm(proposedScomps, mean = MC_chain[i, 5+time+(1:time)], sigma = zigmaS, log = TRUE))
    proposalcurrentScomps<- sum(dmvnorm(MC_chain[i, 5+time+(1:time)], mean = proposedScomps, sigma = zigmaS, log = TRUE))
    
    likelihoodcurrent<- GeneralLoglikelihood(y,MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:time)], MC_chain[i, 5+time+time+(1:ndept)], G(MC_chain[i, 1], MC_chain[i,2]), e_it, MC_chain[i, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)
    likelihoodproposed<- GeneralLoglikelihood(y,MC_chain[i, 5+(1:time)], proposedScomps, MC_chain[i, 5+time+time+(1:ndept)],G(MC_chain[i, 1], MC_chain[i, 2]), e_it, MC_chain[i, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)
    
    mh.ratioS<- exp(likelihoodproposed + priorproposedScomps + proposalcurrentScomps
                    - likelihoodcurrent - priorcurrentScomps - proposalproposedScomps)
    
    print(paste("mh.ratioS =", mh.ratioS))
    
    if(!is.na(mh.ratioS) && runif(1) < mh.ratioS){
      MC_chain[i, 5+time+(1:time)]<- proposedScomps
      
      acceptedS<- acceptedS + 1
    }
    else{
      MC_chain[i, 5+time+(1:time)]<- MC_chain[i, 5+time+(1:time)]
      
      acceptedS<- acceptedS + 0
    }
    
    XnbarPrevR <- XnbarR
    XnbarR <- (i*XnbarR + MC_chain[i, 5+(1:time)])/(i+1)
    zigmaR <- ((i-1)*zigmaR + tcrossprod(MC_chain[i, 5+(1:time)]) + i*tcrossprod(XnbarPrevR) - (i+1)*tcrossprod(XnbarR) + epsilonR*diag(rep(1,time)))/i
    #Robbins Munro tuning
    lambdaR<- lambdaR * exp((2/max(1, i-10)) * (min(mh.ratioR, 1) - 0.234))
    zigmaR<- lambdaR* optconstantR * zigmaR
    #print(zigmaR)
  }   
  
  #Adapting zigmaS
  if(i==10){
    epsilonS<- 0.000007
    XnS<- MC_chain[1:i, 5+time+(1:time)]
    XnbarS <- colMeans(XnS) 
    zigmaS <- cov(XnS) + epsilonS*diag(rep(1, time))
    zigmaS<- optconstantS * zigmaS
  } else if (i > 10){ 
    XnbarPrevS <- XnbarS
    XnbarS <- (i*XnbarS + MC_chain[i, 5+time+(1:time)])/(i+1)
    zigmaS <- ((i-1)*zigmaS + tcrossprod(MC_chain[i, 5+time+(1:time)]) + i*tcrossprod(XnbarPrevS) - (i+1)*tcrossprod(XnbarS) + epsilonS*diag(rep(1,time)))/i
    #Robbins Munro tuning
    lambdaS<- lambdaS * exp((2/max(1, i-10)) * (min(mh.ratioS, 1) - 0.234))
    zigmaS<- lambdaS* optconstantS * zigmaS
    #print(zigmaS)
  }   
  
  #Adapting zigmaU
  if(i==10){
    epsilonU<- 0.000007
    XnU<- MC_chain[1:i, 5+time+time+(1:(ndept-1))]
    XnbarU <- colMeans(XnU) 
    zigmaU <- cov(XnU) + epsilonU*diag(rep(1, ndept-1))
    zigmaU<- optconstantU * zigmaU
  } else if (i > 10){ 
    XnbarPrevU <- XnbarU
    XnbarU <- (i*XnbarU + MC_chain[i, 5+time+time+(1:(ndept-1))])/(i+1)
    zigmaU <- ((i-1)*zigmaU + tcrossprod(MC_chain[i, 5+time+time+(1:(ndept-1))]) + i*tcrossprod(XnbarPrevU) - (i+1)*tcrossprod(XnbarU) + epsilonU*diag(rep(1,ndept-1)))/i
    #Robbins Munro tuning
    lambdaU<- lambdaU * exp((2/max(1, i-10)) * (min(mh.ratioU, 1) - 0.234))
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

#acceptedR/num_iteration
#acceptedS/num_iteration
#accepted/num_iteration
