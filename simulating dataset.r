#This is the code to simulate an example of real data. The dataset has 100 individuals and 10 #SNPs. Data contains binary response variable Y, Sex indicator S (0 for female and 1 for male) #and genotypes of each SNP. Considering minor allele D as the risk allele, dd, dD, DD, d and D #are coded as 0, 1, 2, 0, 2. In this example, No SNP is associated with the binary response Y.
set.seed(1234)
#number of individuals
n=100
#number of SNPs
K=10
#Genotypes
G=matrix(nrow=n,ncol=K)
#Phenotypes
Y=rbinom(n,1,0.5)
#Sex
S=rbinom(n,1,0.5)
#simulate random minor allele frequency between 0 and 0.5
p=runif(K,0,0.5)
for(j in 1:K)
{
   for(i in 1:n)
   {
     if(S[i]==0)
     {
       G[i,j]=sample(c(0,1,2),1,prob=c((1-p[j])^2,2*p[j]*(1-p[j]),p[j]^2))
     }
     else
    {
      G[i,j]=2*rbinom(1,1,p[j])
    }
  }
}
data=cbind(Y,S,G)
write.table(data, "simulated data.txt")