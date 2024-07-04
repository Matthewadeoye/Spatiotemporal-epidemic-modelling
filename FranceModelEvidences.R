source("GeneralLoglikelihood.R")
source("DesignMatrix.R")
source("OutbreakProbability.R")
source("MarginalLikelihood.R")
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


#Model Evidences

#Model 0
fit<- gpu0
draws_df <- fit$draws(format = "df")
posteriordraws<- as.data.frame(draws_df)
fittedKappaU<- posteriordraws[ ,4]
fittedKappaR<- posteriordraws[ ,5]
fittedKappaS<- posteriordraws[ ,6]
fittedUs<- posteriordraws[ ,413:506]
fittedRs<- posteriordraws[ ,100:255]
fittedSs<- posteriordraws[ ,256:411]

PostSamp<- data.frame(KR = fittedKappaR, KS = fittedKappaS, KU = fittedKappaU, 
                      Rs = fittedRs, Ss = fittedSs, Us = fittedUs)
z_it<- DesignMatrixModel0(y, France_adjmat)[[1]]
z_it2<- DesignMatrixModel0(y, France_adjmat)[[2]]
ModelEvidence(y, e_it, R, Model=0, z_it, z_it2, PostSamp, num_samples = 10000)
#Result = -8164.763

#Model 1
fit<- gpu1
draws_df <- fit$draws(format = "df")
posteriordraws<- as.data.frame(draws_df)
fittedG12<- posteriordraws[ ,2]
fittedG21<- posteriordraws[ ,3]
fittedKappaU<- posteriordraws[ ,4]
fittedKappaR<- posteriordraws[ ,5]
fittedKappaS<- posteriordraws[ ,6]
fittedUs<- posteriordraws[ ,414:507]
fittedRs<- posteriordraws[ ,100:255]
fittedSs<- posteriordraws[ ,256:411]
fittedBs<- posteriordraws[ ,412]

PostSamp<- data.frame(G12 = fittedG12, G21 = fittedG21, KR = fittedKappaR, KS = fittedKappaS, 
                      KU = fittedKappaU, Rs = fittedRs, Ss = fittedSs, Us = fittedUs, Bs = fittedBs)
z_it<- DesignMatrixModel1(y, France_adjmat)[[1]]
z_it2<- DesignMatrixModel1(y, France_adjmat)[[2]]
ModelEvidence(y, e_it, R, Model=1, z_it, z_it2, PostSamp, num_samples = 10000)
#Result = -8035.833

#Model 2
fit<- gpu2
draws_df <- fit$draws(format = "df")
posteriordraws<- as.data.frame(draws_df)
fittedG12<- posteriordraws[ ,2]
fittedG21<- posteriordraws[ ,3]
fittedKappaU<- posteriordraws[ ,4]
fittedKappaR<- posteriordraws[ ,5]
fittedKappaS<- posteriordraws[ ,6]
fittedUs<- posteriordraws[ ,414:507]
fittedRs<- posteriordraws[ ,100:255]
fittedSs<- posteriordraws[ ,256:411]
fittedBs<- posteriordraws[ ,412]

PostSamp<- data.frame(G12 = fittedG12, G21 = fittedG21, KR = fittedKappaR, KS = fittedKappaS, 
                      KU = fittedKappaU, Rs = fittedRs, Ss = fittedSs, Us = fittedUs, Bs = fittedBs)
z_it<- DesignMatrixModel2(y, France_adjmat)[[1]]
z_it2<- DesignMatrixModel2(y, France_adjmat)[[2]]
ModelEvidence(y, e_it, R, Model=2, z_it, z_it2, PostSamp, num_samples = 10000)
#Result = -8058.945

#Model 3
fit<- gpu3
draws_df <- fit$draws(format = "df")
posteriordraws<- as.data.frame(draws_df)
fittedG12<- posteriordraws[ ,2]
fittedG21<- posteriordraws[ ,3]
fittedKappaU<- posteriordraws[ ,4]
fittedKappaR<- posteriordraws[ ,5]
fittedKappaS<- posteriordraws[ ,6]
fittedUs<- posteriordraws[ ,415:508]
fittedRs<- posteriordraws[ ,100:255]
fittedSs<- posteriordraws[ ,256:411]
fittedBs<- posteriordraws[ ,412:413]

PostSamp<- data.frame(G12 = fittedG12, G21 = fittedG21, KR = fittedKappaR, KS = fittedKappaS, 
                      KU = fittedKappaU, Rs = fittedRs, Ss = fittedSs, Us = fittedUs, Bs = fittedBs)
z_it<- DesignMatrixModel3(y, France_adjmat)[[1]]
z_it2<- DesignMatrixModel3(y, France_adjmat)[[2]]
ModelEvidence(y, e_it, R, Model=3, z_it, z_it2, PostSamp, num_samples = 10000)
#Result = -8061.628

#Model 4
fit<- gpu4
draws_df <- fit$draws(format = "df")
posteriordraws<- as.data.frame(draws_df)
fittedG12<- posteriordraws[ ,2]
fittedG21<- posteriordraws[ ,3]
fittedKappaU<- posteriordraws[ ,4]
fittedKappaR<- posteriordraws[ ,5]
fittedKappaS<- posteriordraws[ ,6]
fittedUs<- posteriordraws[ ,414:507]
fittedRs<- posteriordraws[ ,100:255]
fittedSs<- posteriordraws[ ,256:411]
fittedBs<- posteriordraws[ ,412]

PostSamp<- data.frame(G12 = fittedG12, G21 = fittedG21, KR = fittedKappaR, KS = fittedKappaS, 
                      KU = fittedKappaU, Rs = fittedRs, Ss = fittedSs, Us = fittedUs, Bs = fittedBs)
z_it<- DesignMatrixModel4(y, France_adjmat)[[1]]
z_it2<- DesignMatrixModel4(y, France_adjmat)[[2]]
ModelEvidence(y, e_it, R, Model=4, z_it, z_it2, PostSamp, num_samples = 10000)
#Result = -8080.507

#Model 5
fit<- gpu5
draws_df <- fit$draws(format = "df")
posteriordraws<- as.data.frame(draws_df)
fittedG12<- posteriordraws[ ,2]
fittedG21<- posteriordraws[ ,3]
fittedKappaU<- posteriordraws[ ,4]
fittedKappaR<- posteriordraws[ ,5]
fittedKappaS<- posteriordraws[ ,6]
fittedUs<- posteriordraws[ ,414:507]
fittedRs<- posteriordraws[ ,100:255]
fittedSs<- posteriordraws[ ,256:411]
fittedBs<- posteriordraws[ ,412]

PostSamp<- data.frame(G12 = fittedG12, G21 = fittedG21, KR = fittedKappaR, KS = fittedKappaS, 
                      KU = fittedKappaU, Rs = fittedRs, Ss = fittedSs, Us = fittedUs, Bs = fittedBs)
z_it<- DesignMatrixModel5(y, France_adjmat)[[1]]
z_it2<- DesignMatrixModel5(y, France_adjmat)[[2]]
ModelEvidence(y, e_it, R, Model=5, z_it, z_it2, PostSamp, num_samples = 10000)
#Result = -8076.989

#Model 6
fit<- gpu6
draws_df <- fit$draws(format = "df")
posteriordraws<- as.data.frame(draws_df)
fittedG12<- posteriordraws[ ,2]
fittedG21<- posteriordraws[ ,3]
fittedKappaU<- posteriordraws[ ,4]
fittedKappaR<- posteriordraws[ ,5]
fittedKappaS<- posteriordraws[ ,6]
fittedUs<- posteriordraws[ ,415:508]
fittedRs<- posteriordraws[ ,100:255]
fittedSs<- posteriordraws[ ,256:411]
fittedBs<- posteriordraws[ ,412:413]

PostSamp<- data.frame(G12 = fittedG12, G21 = fittedG21, KR = fittedKappaR, KS = fittedKappaS, 
                      KU = fittedKappaU, Rs = fittedRs, Ss = fittedSs, Us = fittedUs, Bs = fittedBs)
z_it<- DesignMatrixModel6(y, France_adjmat)[[1]]
z_it2<- DesignMatrixModel6(y, France_adjmat)[[2]]
ModelEvidence(y, e_it, R, Model=6, z_it, z_it2, PostSamp, num_samples = 10000)
#Result = -8111.527