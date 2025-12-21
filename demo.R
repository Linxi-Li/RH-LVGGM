# Load required functions
# install.packages(c("MASS", "glasso", "flexclust"))
source("main.R")
source("sim.R")
source("comparison.R")

#  Data Generation 
data.gener=Generate.data(n=1500,p=40,r=5,pi0=0.04,piK=c(0.32,0.32,0.32),K=3,sparsity=0.1,mue=1.5,sparsity.latent=0.7,seed=1) 
data=data.gener$data

# Generate test data (same parameters and seed)
data.test.gener=Generate.data(n=500,p=40,r=5,pi0=0.04,piK=c(0.32,0.32,0.32),K=3,sparsity=0.1,mue=1.5,sparsity.latent=0.7,seed=1) 
data.test=data.test.gener$data

# Define regularization parameter search space
lambda1.seq=seq(0.001,0.1,length.out =15)
lambda2.seq=seq(0.1,1,length.out =15)


# ===== ===== ===== Proposed Method ===== ===== =====
# Initialize clustering with noise points
tau.ini=InitClust.nosie(data,G=3)

# Fitting model
fit.our=OT.RH.lvggm(data,G=3,tau.ini,lambda1.seq,lambda2.seq,npr.max=0.5, erc=50, iter.max=100, tol=1e-03)

# Performance evaluation
Metric(data.test,data.test.gener$member,fit.our$lambda.1,fit.our$lambda.2,fit.our,data.gener$S,data.gener$Theta,fit.our$S,fit.our$Theta,fit.our$L,data.gener$L,data.gener$member,fit.our$cluster)


# ===== ===== ===== Comparison Methods ===== ===== =====

#  LVGGM
fit.lvggm=Compare.lvggm(data,lambda1.seq,lambda2.seq)

# H-LVGGM
tau.ini=InitClust(data,G=3)
fit.h.lvggm=Compare.H.lvggm(data,G=3,tau.ini,lambda1.seq,lambda2.seq,npr.max=0.5, erc=50, iter.max=100, tol=1e-03)

# OR-LVGGM
fit.or.lvggm=Compare.or.lvggm(data,K=3,data.gener$member,lambda1.seq,lambda2.seq)

#  Glasso
fit.glasso=Compare.glasso(data,lambda1.seq)

# OR-Glasso
fit.orglasso=Compare.or.glasso(data,K=3,data.gener$member,lambda1.seq)

#  RH-GGM
tau.ini=InitClust.nosie(data,G=3)
fit.rh.ggm=OT.RH.lvggm(data,G=3,tau.ini,lambda1.seq,lambda2.seq=1000,npr.max=0.5, erc=50, iter.max=100, tol=1e-03)

#  H-GGM
tau.ini=InitClust(data,G=3)
fit.h.ggm=Compare.H.lvggm(data,G=3,tau.ini,lambda1.seq,lambda2.seq=1000,npr.max=0.5, erc=50, iter.max=100, tol=1e-03)