#This code performs frequentist and BMA analysis on the simualted dataset
library(MCMCpack)
library(coda)
library(mvtnorm)
data=read.table("simulated data.txt",header=T)

#number of MCMC samples
M=1000

#precision parameter lambda
precision=1

#number of individuals
n=dim(data)[1]

#number of SNPs
K=dim(data)[2]-2

MAF=rep(0,K)

#effect sizes and p-values under XCI model and no XCI model
beta1=rep(0,K)
beta2=rep(0,K)
P1=rep(0,K)
P2=rep(0,K)

#p-value of Z_max
min.P=rep(0,K)

BF12=rep(0,K)
BFAN=rep(0,K)
Lower=rep(0,K)
Upper=rep(0,K)

for(k in 1:K)
{
  Y=data[,1]
  S=data[,2]
  G=data[,k+2]

#Define genotype coding under XCI and no XCI assumptions for each SNP
  G1=G/2
  G2=G-S
  G2[G2==-1]=0

#computing minor allele frequencies
  male=G[S==1]
  female=G[S==0]
  MAF[k]=(length(male[male==2])+length(female[female==2])+length(female[female==1])/2)/length(S)

#frequentist inferences of beta and p-value under XCI and no XCI model, and p-value of Z_max
  fit1=glm(Y~G1+S,family="binomial")
  fit2=glm(Y~G2+S,family="binomial")
  fit0=glm(Y~S,family="binomial")
  beta1[k]=summary(fit1)$coefficients[2,4]
  beta2[k]=summary(fit2)$coefficients[2,4]
  Z1=summary(fit1)$coefficients[2,3]
  Z2=summary(fit2)$coefficients[2,3]
  P1[k]=summary(fit1)$coefficients[2,4]
  P2[k]=summary(fit2)$coefficients[2,4]

#theoretically adjusted p-value for Z_max
  z=max(abs(Z1),abs(Z2))
  r=cor(G1,G2)
  min.P[k]=1-pmvnorm(lower=c(-z,-z), upper=c(z,z),  sigma = cbind(c(1,r),c(r,1)))

#Baeysian inference
  BF12[k]=exp((BIC(fit2)-BIC(fit1))/2)
  BF1N=exp((BIC(fit0)-BIC(fit1))/2)
  BF2N=exp((BIC(fit0)-BIC(fit2))/2)
  BFAN[k]=(BF1N+BF2N)/2
  X1=cbind(rep(1,n),G1,S)
  X2=cbind(rep(1,n),G2,S)

#simulating M posterior samples under each model
  pos1=MCMClogit(Y~G1+S,B0=t(X1)%*%X1*precision/n,mcmc=M)
  pos2=MCMClogit(Y~G2+S,B0=t(X2)%*%X2*precision/n,mcmc=M)

#draw posterior samples of BMA model
  pos=matrix(nrow=M,ncol=3)
  for(j in 1:M)
  {
    I=rbinom(1,1,BF12[k]/(1+BF12[k]))
    if(I==1)
    {
      pos[j,]=pos1[j,]
    }
    else
    {
      pos[j,]=pos2[j,]
    }
  }

#compute the BMA-based HPD interval 
  mc=as.mcmc(pos[,2])
  MC.HPD=HPDinterval(mc)
  Lower[k]=MC.HPD[1]
  Upper[k]=MC.HPD[2]
}
#Summary result analogous to Table 1 in the paper
result=cbind(MAF,beta1,beta2,P1,P2,min.P,Lower,Upper,BF12,BFAN)

#Summary plot analogous to Figure 3 in the paper. SNPs ranked by the lower bound of HPD
#intervals
jpeg('SNPs.jpeg',width=1500,height=600)
par(mfrow=c(1,3),cex=1, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, mai = c(0.9, 0, 0.5, 0), omi = c(0, 1, 0, 0.1))
y=c(1:K)
order=rev(order(Lower))
plot(NULL,ylim=rev(range(c(1,K))),xlim=range(c(0,-log10(min.P))), xlab="-log10 p-values of Z_max", axes = FALSE)
axis(1)
axis(2)
box()
arrows(0, y, -log10(min.P[order]), y, length=0.01, angle=90, code=3)

plot(NULL,ylim=rev(range(c(1,K))),xlim=range(c(0,log10(BFAN))), xlab="log10 BF_AN", axes = FALSE)
axis(1)
box()
arrows(0, y, log10(BFAN[order]), y, length=0.01, angle=90, code=3)

plot(NULL,ylim=rev(range(c(1,K))),xlim=range(c(Lower, Upper)), xlab="BMA HPD intervals", axes=FALSE)
 axis(1)
 box()
arrows(Lower[order], y, Upper[order], y, length=0.01, angle=90, code=3)

mtext("SNPs", side = 2, outer = TRUE, cex = 2, line = 3)
dev.off()
