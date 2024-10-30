
#install.packages("kmeansstep", repos="http://R-Forge.R-project.org")

library(simstudy)
library(mclust)
library(cluster)
library(fpc)
library(clValid)
library(kmeansstep)
library(fossil)
library(clusterSim)
library(cstab)


mupe = function(p,df) exp(pchisq(qchisq(p,df=df),df=df+2,log.p=T) - log(p))*df




dabic = function(x,cluster){
	p = ncol(x)
	n = nrow(x)
	pen = 2*p*log(n)

	lvup = log(diag(var(x)*(n-1)/n))
	nlik = 0.5*n*sum(lvup)
	bic = bica = nlik + 0.5*( pen )

	nn = table(cluster)
	kk = length(nn)
	if(kk>1){
		pr = nn/n
		mp = pr
		for(j in 1:length(mp)) mp[j] = mupe(pr[j],df=1)
		msq = msqc = numeric(length(pr))
		for(jj in 1:kk){
			if(nn[jj]==1){
				msq[jj] = sum(lvup)
				msqc[jj] = sum(lvup)
			}else{
				msqjj = log(diag(var(x[cluster==jj,,drop=FALSE])*(nn[jj]-1)/nn[jj]))
				msqjjc = msqjj - log(mp[jj])
				if(any(msqjjc>lvup)) msqjjc[msqjjc>lvup] = lvup[msqjjc>lvup]
				msq[jj] = sum(msqjj)
				msqc[jj] = sum(msqjjc)
			}
		}
		nLik = 0.5*( sum(msq*pr)*n )
		Bic = nLik + 0.5*( pen*kk )
		Bica = 0.5*( sum(msqc*pr)*n + pen*kk )
	}else{
		nLik = nlik
		Bic = bic
		Bica = bica
	}

	return( list(nlik=nLik,bic=Bic,bica=Bica) )
}




selectKmeans = function(x,Kmax=10,Kplot=2){
	p = ncol(x); n = nrow(x)
	pen = 2*p*log(n)
	bic = bica = bicn = sil = chi = numeric(Kmax)
	clusterbic = clusterbica = clusterbicn = clustersil = clusterchi = rep(1,n)

	da1 = dabic(x,rep(1,n))
	bic[1] = da1$bic
	bica[1] = da1$bica	
	sil[1] = chi[1] = -Inf
	bicn[1] = kmeansBIC(kmeans(scale(x),centers=1,nstart=25))

	for(kk in 2:Kmax){
		km = kmeans(scale(x),centers=kk,nstart=25)
		dakk = dabic(x,km$cluster)
		if(kk==Kplot) plot(x,col=km$cluster)
		bic[kk] = dakk$bic
		bica[kk] = dakk$bica
		ss = silhouette(km$cluster, dist(scale(x)))
		sil[kk] = mean(ss[, 3])
		bicn[kk] = kmeansBIC(km)
		chi[kk] = calinhara(x,km$cluster)
		if(min(bic[1:(kk-1)])>bic[kk]) clusterbic = km$cluster
		if(min(bica[1:(kk-1)])>bica[kk]) clusterbica = km$cluster
		if(max(sil[1:(kk-1)])<sil[kk]) clustersil = km$cluster
		if(max(chi[1:(kk-1)])<chi[kk]) clusterchi = km$cluster
		if(min(bicn[1:(kk-1)])>bicn[kk]) clusterbicn = km$cluster
	}
	return(list(bic=bic,bica=bica,sil=sil,chi=chi,bicn=bicn,clusterbic=clusterbic,clusterbica=clusterbica,clusterbicn=clusterbicn,clustersil=clustersil,clusterchi=clusterchi))
}






