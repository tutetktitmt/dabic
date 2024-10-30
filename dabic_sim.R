

library(simstudy)
library(mclust)
library(cluster)
library(fpc)
library(clValid)
library(kmeansstep)
library(aricode)
library(clusterSim)
library(cstab)



source("dabic.R")



################################################

set.seed(71042)


K = 1  # c(1,2,3,5,10)
n = 1000  # c(1000,2000)
rho = -0.0  # c(-0.0,-0.7)
p = 2  # c(2,3,5)
dist = "normal"  # c("normal","uniform")
vconstK = T  # c(T,F)

Kmax = 20
R = 50

args = commandArgs(trailingOnly=TRUE)

print(length(args))
print(args)

if(length(args)==7){
	# Para = expand.grid( c(1,2,3,5,10), c(1000,2000), c(2,3,5))
	# Para[l,1]," ",Para[l,2]," ",Para[l,3]," ",dist," ",rho," ",vconstK," ",R
	K = as.numeric(args[1])
	n = as.numeric(args[2])
	p = as.numeric(args[3])
	dist = args[4]
	rho = as.numeric(args[5])
	vconstK = as.numeric(args[6])
	R = as.numeric(args[7])
	# "3"      "1000"   "2"      "normal" "-0.7"   "0"      "500"
	fi = paste0("simresult/cl-abic",K,"_",n,"_",p,"_",dist,"_",rho,"_",vconstK,"_",R,".RData")
	print(fi)
}


xmu00 = matrix(0,p,K)
sig200 = K:1
if(vconstK) sig200[] = 2
if(dist=="normal") sig200 = sig200/8
xSig00 = rho^abs(outer(1:p,1:p,"-"))

xmu00[,1] = 3
if(K>1) for(i in 2:K) xmu00[,i] = 3*i


################################################

OptK = matrix(0,R,9); colnames(OptK) = c("daBIC","BICcomv1","naBIC","CHI","Silhouette","GapStat","JumpStat","GMBIC","PS")
ARand = OptK

for(r in 1:R){

#
im = sample(1:K,n,replace=T)

x = genCorGen(n, nvars=p, params1=xmu00[,1], params2=sig200[1]+(dist=="uniform")*xmu00[,1], dist=dist, corMatrix=xSig00, wide=T)[,-1]

if(K>1){
	for(i in 2:K) x[im==i,] = genCorGen(sum(im==i), nvars=p, params1=xmu00[,i], params2=sig200[i]+(dist=="uniform")*xmu00[,i], dist=dist, corMatrix=xSig00, wide=T)[,-1]
}

plot(x)
#




# adjusted BIC
abica = selectKmeans(x,Kmax=Kmax,Kplot=K)
OptK[r,"naBIC"] = which(min(abica$bic,na.rm=T)==abica$bic)
OptK[r,"daBIC"] = which(min(abica$bica,na.rm=T)==abica$bica)
ARand[r,"naBIC"] = ARI(im,abica$clusterbic)
ARand[r,"daBIC"] = ARI(im,abica$clusterbica)


# Silhouette score
OptK[r,"Silhouette"] = which(max(abica$sil,na.rm=T)==abica$sil)
ARand[r,"Silhouette"] = ARI(im,abica$clustersil)


#  Calinski-Harabasz index
OptK[r,"CHI"] = which(max(abica$chi,na.rm=T)==abica$chi)
ARand[r,"CHI"] = ARI(im,abica$clusterchi)


# kmeans BIC
OptK[r,"BICcomv1"] = which(min(abica$bicn,na.rm=T)==abica$bicn)
ARand[r,"BICcomv1"] = ARI(im,abica$clusterbicn)


# Gap stat
resg = clusGap(scale(x),kmeans,K.max=Kmax,B=100)
optresg = with(resg,maxSE(Tab[,"gap"],Tab[,"SE.sim"]))
OptK[r,"GapStat"] = optresg
kmresg = kmeans(scale(x),optresg,nstart=25)
ARand[r,"GapStat"] = ARI(im,kmresg$cluster)


# Sugar-James Jump stat
cDst = cDistance(scale(x), kseq=2:Kmax)
OptK[r,"JumpStat"] = cDst$k_Jump
kmresj = kmeans(scale(x),cDst$k_Jump,nstart=25)
ARand[r,"JumpStat"] = ARI(im,kmresj$cluster)


# Gauss mixture BIC
#mclbic = mclustBIC(scale(x),G=1:Kmax)
mclbic = mclustBIC(x,G=1:Kmax)
summary(mclbic)
OptK[r,"GMBIC"] = which(mclbic==max(mclbic,na.rm=T),arr.ind=T)[1,1]
#mclopt = mclustModel(scale(x), mclbic)
mclopt = mclustModel(x, mclbic)
ARand[r,"GMBIC"] = ARI(im,apply(mclopt$z,1,which.max))


# prediction strength
prst = prediction.strength(scale(x),2,Kmax,M=3)
OptK[r,"PS"] = prst$optimalk
kmprst = kmeans(scale(x),prst$optimalk,nstart=25)
ARand[r,"PS"] = ARI(im,kmprst$cluster)


}


OptK
ARand



save.image(fi)




