#differential correlaion test using Fisher's Z transformation
#http://www.sciencedirect.com/science/article/pii/S0378111912014497
dcor=function(r1,r2,n1,n2,two.tailed=FALSE){
	if(any(abs(r1)>1) || any(abs(r2)>1)) stop('Invalid correlation value(s)\n')
	if(any(n1 <= 3) || any(n2 <= 3)) stop('n1 and n2 must be larger than 3\n')
	stopifnot(length(r1) == length(r2))
	r1[r1 > 0.999999999] = 0.999999999 #this is set to avoid log(0)
	r1[r1 < -0.999999999] = -0.999999999
	r2[r2 > 0.999999999] = 0.999999999
	r2[r2 < -0.999999999] = -0.999999999
	z1=log((1+r1)/(1-r1))/2
	z2=log((1+r2)/(1-r2))/2
	z=(z1-z2)/sqrt(1/(n1-3)+1/(n2-3))
	p=pnorm(abs(z),lower.tail = FALSE)
	if(two.tailed) p=2*p
	if(length(z)>1) return(list(z=z,p=p))
	c(z=z,p=p)
}
r2p=function(r,n,alternative=c('two.sided','less','greater')){
	alternative=match.arg(alternative)
	df <- n - 2L
	STATISTIC <- c(t = sqrt(df) * r/sqrt(1 - r^2)) #how to handle the case when r^2 equals to 1?
	switch(alternative, less = pt(STATISTIC, df), 
	greater = pt(STATISTIC, df, lower.tail = FALSE), 
	two.sided = 2 * pmin(pt(STATISTIC, df), pt(STATISTIC, df, lower.tail = FALSE)))
}
setGeneric(name="ddcor",
	def=function(cor1, cor2, n1, n2, two.tailed=FALSE, adjust='BH', FDR=0.05,...){
		standardGeneric("ddcor")
	},signature=c("cor1", "cor2")
)
setMethod(f = "ddcor",
	signature=c(cor1 = "matrix", cor2 = "matrix"),
	definition = function(cor1, cor2, n1, n2, two.tailed=FALSE, adjust='BH', FDR=0.05,selection=c('whole','upper.tri','lower.tri'),diag=FALSE, ...){
	selection=match.arg(selection)
	stopifnot(identical(dim(cor1),dim(cor2)))
	if(is.null(rownames(cor1)) || is.null(colnames(cor1))) stop('cor1 must have rownames and colum names\n')
	stopifnot(identical(dimnames(cor1),dimnames(cor2)))
	#
	if(selection=='upper.tri'){
		p=dcor(cor1[upper.tri(cor1)], cor2[upper.tri(cor2)], n1, n2, two.tailed=two.tailed)
	}else if(selection=='lower.tri'){
		p=dcor(cor1[lower.tri(cor1)], cor2[lower.tri(cor2)], n1, n2, two.tailed=two.tailed)
	}else{
		p=dcor(cor1, cor2, n1, n2, two.tailed=two.tailed)
	}
	z=as.vector(p$z)
	p=as.vector(p$p)
	Q=p.adjust(p,method=adjust)
	p1=as.vector(r2p(cor1,n1))
	p2=as.vector(r2p(cor2,n2))
	q1=p.adjust(p1,method=adjust)
	q2=p.adjust(p2,method=adjust)
	#
	sig = which(Q<=FDR)
	if(length(sig)==0) stop('No significance found at given FDR threshold.\n')
	#
	if(selection=='upper.tri'){
		ddres = df0(Var=arrayInd(which(upper.tri(cor1))[sig],.dim=dim(cor1)),Cor1=cor1[upper.tri(cor1)][sig],Cor2=cor2[upper.tri(cor2)][sig])
	}else if(selection=='lower.tri'){
		ddres = df0(Var=arrayInd(which(lower.tri(cor1))[sig],.dim=dim(cor1)),Cor1=cor1[lower.tri(cor1)][sig],Cor2=cor2[lower.tri(cor2)][sig])
	}else{
		ddres = df0(Var=arrayInd(sig,.dim=dim(cor1)),Cor1=as.vector(cor1[sig]),Cor2=as.vector(cor2[sig]))
	}
	#
	ddres = df0(ddres,
		p1=p1[sig],p2=p2[sig],
		q1=q1[sig],q2=q2[sig],
		z=z[sig],p=p[sig],p.adj=Q[sig]
	)
	ddres$Classes=paste(ifelse(ddres$q1 <= FDR,ifelse(ddres$Cor1 < 0,'-','+'),'0'),ifelse(ddres$q2 <= FDR,ifelse(ddres$Cor2 < 0,'-','+'),'0'),sep='/')
	ddres$Var.1=rownames(cor1)[ddres$Var.1]
	ddres$Var.2=colnames(cor1)[ddres$Var.2]
	colnames(ddres)[1:2]=c('Var1','Var2')
	ddres
})
setMethod(f = "ddcor",
	signature=c(cor1 = "numeric", cor2 = "numeric"),
	definition = function(cor1, cor2, n1, n2, two.tailed=FALSE, adjust='BH', FDR=0.05,...){
	stopifnot(length(cor1) == length(cor2))
	stopifnot(identical(names(cor1),names(cor2)))
	#
	p=dcor(cor1, cor2, n1, n2, two.tailed=two.tailed)
	z=as.vector(p$z)
	p=as.vector(p$p)
	Q=p.adjust(p,method=adjust)
	p1=as.vector(r2p(cor1,n1))
	p2=as.vector(r2p(cor2,n2))
	q1=p.adjust(p1,method=adjust)
	q2=p.adjust(p2,method=adjust)
	#
	sig = which(Q<=FDR)
	if(length(sig)==0) stop('No significance found at given FDR threshold.\n')
	ddres = df0(Var=sig,
		Cor1=cor1[sig],Cor2=cor2[sig],
		p1=p1[sig],p2=p2[sig],
		q1=q1[sig],q2=q2[sig],
		z=z[sig],p=p[sig],p.adj=Q[sig]
	)
	ddres$Classes=paste(ifelse(ddres$q1 <= FDR,ifelse(ddres$Cor1 < 0,'-','+'),'0'),ifelse(ddres$q2 <= FDR,ifelse(ddres$Cor2 < 0,'-','+'),'0'),sep='/')
	if(!is.null(names(cor1))) ddres$Var=names(cor1)[ddres$Var]
	ddres
})
MDC=function(MatA,MatB,genes,modules,cor.method = c("pearson", "spearman"),dCorAvgMethod=c('mean','median','wilcox'),cor2connectivity=ifelse(dCorAvgMethod=='wilcox',TRUE,FALSE),intType=ifelse(cor2connectivity,1,0),Power=1,adjust.method='BH',nPerms=1000,random.seed=12345,ncores=1){
	cor.method=match.arg(cor.method)
	dCorAvgMethod=match.arg(dCorAvgMethod)
	stopifnot(nrow(MatA)==nrow(MatB))
	stopifnot(identical(rownames(MatA),rownames(MatB)))
	stopifnot(all(genes %in% rownames(MatA)))
	if(ncol(MatA)<3 || ncol(MatB) <3) stop('Too few samples in the input data\n')
	if(any(is.na(MatA)) || any(is.na(MatB))) stop('Missing data not allowed in the input data\n')
	MatA=MatA[rownames(MatA) %in% genes,]
	MatB=MatB[rownames(MatA),]
	nA=ncol(MatA)
	nB=ncol(MatB)
	uniqueLabels=unique(modules)
	cl=NULL
	if(ncores>1 && require('parallel')){
		cl=makeCluster(ncores)
		cat('Use',ncores,'cores\n')
	}
	diffcors=.MDCcaller1(MatA=MatA,MatB=MatB,uniqueLabels=uniqueLabels,modules=modules,genes=genes,cor.method=cor.method,dCorAvgMethod=dCorAvgMethod,cor2connectivity=cor2connectivity,intType=intType,Power=Power,diff.only=FALSE,cl=cl)
	if(dCorAvgMethod=='wilcox'){
		padj=as.vector(p.adjust(diffcors[,'P.Value'],method=adjust.method))
		res=data.frame(Module=uniqueLabels,diffcors,adj.P.Value=padj,stringsAsFactors=FALSE)
	}else{
		cat('Start permutation analysis.\n')
		MatC=cbind(MatA,MatB)
		pvals=rep(NA,length(uniqueLabels))
		pvals[]=0
		set.seed(random.seed)
		for(j in 1:nPerms){
			if( j %% 10 ==0) cat('Running permutation',j,'...\n')
			s1=sample(ncol(MatC),nA)
			dds=.MDCcaller1(MatA=MatC[,s1],MatB=MatC[,-s1],uniqueLabels=uniqueLabels,modules=modules,genes=genes,cor.method=cor.method,dCorAvgMethod=dCorAvgMethod,cor2connectivity=cor2connectivity,intType=intType,Power=Power,diff.only=TRUE,cl=cl)
			pvals = pvals + as.numeric(abs(dds) >= abs(diffcors[,'Diff']) )
		}
		pvals=(pvals+1)/(nPerms+1)
		padj=p.adjust(pvals,method=adjust.method)
		res=data.frame(Module=uniqueLabels,diffcors,P.Value=pvals,adj.P.Value=padj,stringsAsFactors=FALSE)
	}
	if(!is.null(cl)) stopCluster(cl)
	res
}
.MDCcaller1=function(MatA,MatB,uniqueLabels,modules,genes,cor.method=c("pearson", "spearman"),dCorAvgMethod=c('mean','median','wilcox'),cor2connectivity=ifelse(dCorAvgMethod=='wilcox',TRUE,FALSE),intType=ifelse(cor2connectivity,1,0),Power=1,diff.only=TRUE,cl=NULL){
	cor.method=match.arg(cor.method)
	if(cor.method=='spearman'){
		MatA=t(apply(MatA,1,rank))
		MatB=t(apply(MatB,1,rank))
	}
	MatA=scale(t(MatA))
	MatB=scale(t(MatB))
	func1=function(x,uniqueLabels,modules,genes,MatA,MatB,dCorAvgMethod,cor2connectivity,intType,Power,diff.only){
		g=genes[modules==uniqueLabels[x]]
		if(length(g)==1){
			if(diff.only) return(NA)
			if(dCorAvgMethod=='mean' || dCorAvgMethod=='median'){
				return(c(meanA=NA,meanB=NA,medianA=NA,medianB=NA,Diff=NA))
			}
			return(c(meanA=NA,meanB=NA,medianA=NA,medianB=NA,w=NA,z=NA,P.Value=NA))
		}
		.moduleDiffCon(MatA=MatA[,g],MatB=MatB[,g],intType=intType,Power=Power,cor2connectivity=cor2connectivity,dCorAvgMethod=dCorAvgMethod,diff.only=diff.only)
	}
	if(is.null(cl)){
		diffcors=lapply(1:length(uniqueLabels),func1,uniqueLabels=uniqueLabels,modules=modules,genes=genes,MatA=MatA,MatB=MatB,dCorAvgMethod=dCorAvgMethod,cor2connectivity=cor2connectivity,intType=intType,Power=Power,diff.only=diff.only)
	}else{
		a=clusterExport(cl,varlist=c('.moduleDiffCon','dat2adjacency'),envir=environment())
		diffcors=parLapply(cl,1:length(uniqueLabels),func1,uniqueLabels=uniqueLabels,modules=modules,genes=genes,MatA=MatA,MatB=MatB,dCorAvgMethod=dCorAvgMethod,cor2connectivity=cor2connectivity,intType=intType,Power=Power,diff.only=diff.only)
	}
	if(diff.only==FALSE || dCorAvgMethod=='wilcox'){
		diffcors=do.call(rbind,diffcors)
	}else{
		diffcors=unlist(diffcors)
	}
	diffcors
}
.moduleDiffCon=function(MatA,MatB,dCorAvgMethod=c('mean','median','wilcox'),cor2connectivity=ifelse(dCorAvgMethod=='wilcox',TRUE,FALSE),intType=ifelse(cor2connectivity,1,0),Power=1,diff.only=TRUE){
#MatA and MatB, columns are standardized variables and rows ae samples
	dCorAvgMethod=match.arg(dCorAvgMethod)
	stopifnot(ncol(MatA)==ncol(MatB))
	nvar=ncol(MatA)
	corA=dat2adjacency(MatA,intType=intType,Power=Power)
	corB=dat2adjacency(MatB,intType=intType,Power=Power)
	stopifnot(nrow(corA)==nrow(corB))
	if(cor2connectivity){
		corA=rowSums(corA,na.rm=TRUE)/(nvar-1)
		corB=rowSums(corB,na.rm=TRUE)/(nvar-1)
	}else{
		corA=corA[upper.tri(corA)]
		corB=corB[upper.tri(corB)]
	}
	res=c()
	if(!diff.only){
		res=c(meanA=mean(corA,na.rm=TRUE),meanB=mean(corB,na.rm=TRUE),medianA=median(corA,na.rm=TRUE),medianB=median(corB,na.rm=TRUE))
	}
	if(dCorAvgMethod=='mean') return(c(res,Diff=mean(corA-corB,na.rm=TRUE)))
	if(dCorAvgMethod=='median') return(c(res,Diff=median(corA-corB,na.rm=TRUE)))
	#Wilcox test
	fit=wilcox.test(x=corA,y=corB,paired = TRUE, exact=FALSE)
	wilcox.approx.z=(fit$statistic[['V']]-nvar*(nvar+1)/4)/sqrt(nvar*(nvar+1)*(2*nvar-1)/24)
	c(res,w=fit$statistic[['V']],nvar=nvar,z=wilcox.approx.z,P.Value=fit$p.value)
}
dat2adjacency=function(x,intType=1,Power=1){
	corMat=t(x) %*% x/(nrow(x)-1)
	stopifnot(nrow(corMat)==ncol(x))
	diag(corMat)=NA
	if(intType==1){
		corMat=abs(corMat)
	} else if (intType==2){
		corMat = (1+corMat)/2;
	} else if (intType==3){
		corMat[corMat < 0] = 0; 
	}
	corMat^Power
}
ratio.test=function(x,y,covxy=0){
#reference: A biologist's guide to statistical thinking and analysis" by David S. Fay and Ken Gerow
#http://www.wormbook.org/chapters/www_statisticalanalysis/statisticalanalysis.html
	se=function(x){
		x=deEmptyNA(x,drop.empty.string = FALSE)
		sd(x)/sqrt(length(x))
	}
	mx=mean(x,na.rm=TRUE)
	my=mean(y,na.rm=TRUE)
	ratio=mx/my
	se.ratio=ratio*sqrt((se(x)/mx)^2+(se(y)/my)^2-2*covxy/mx/my)
	c(ratio=ratio,se=se.ratio)
}
#Computes the Mann-Kendall's trend test
MannKendall.test=function(x){
#x is a vector of time series data; ordered by time
#Note that P-value is obtained from normal approximation.
	n=length(x)
	if(n<2) stop('Number of data points is too small\n')
	T = ((n-1)*n) / 2
	S = sum( sapply( 2:n, function(i) sum( sign(x[i]-x[1:(i-1)]) ) ) )
	tau = S / T
	varS = n * (n-1) * (2*n + 5) / 18
	S=S+ifelse(S>0,-1,1) #continuity correction
	p = 2 * pnorm(-abs(S)/sqrt(varS)) 	#two-sided p value
	return(c(tau=tau,p.value=p))
}
#Computes Jonckheere's non-parametric trend test for ordered alternatives
Jonckheere.test=function(y,x){
#Input: y, a vector of numeric values; x, a vector or an ordered factor specifying the grouping of each values in y
#Reference: Hollander and Wolfe. 1999. Nonparametric statistical methods (Second edition). John Wiley & Sons, Inc., New York.
#Output: a vector with three named elements, J (the test statistics), Z (normalized J value) and p.value (two-sided p value calculated from normal approximation)
	is.rm=is.na(y) | is.na(x) # NA will be removed
	y=y[!is.rm]
	x=x[!is.rm]
	dat=split(y,as.ordered(x))
	dat=dat[sapply(dat,length)>0]
	M=length(dat)
	if(M<3) stop('At least 3 groups are required\n')
	nt=sapply(dat,function(a) length(a))
	dat=unlist(dat)
	cnt=c(0,cumsum(nt))
	N=length(y)
	J=sapply(1:(M-1),function(i){
		sum(sapply(dat[(cnt[i]+1):cnt[i+1]],function(d) sum(dat[(cnt[i+1]+1):N]>d)))
	})
	J=sum(J)
	eJ=(N*N-sum(nt*nt))/4
	varJ=(N*N*(2*N+3)-sum(nt*nt*(2*nt+3)))/72
	Z=(J-eJ)/sqrt(varJ)
	p=2*pnorm(-abs(Z))
	return(c(J=J,Z=Z,p.value=p))
}
#Jensen-Shannon divergence
JSD <- function(mat, pseudocount=0.000001) {
	KLD <- function(x,y) sum(x *log(x/y))
	JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
	mat[mat==0]=pseudocount
	result=apply(mat,2,function(x,dat) apply(dat,2,function(y,x) JSD(x,y),x=x),dat=mat)
	colnames(result) = rownames(result) = colnames(mat)
	result=as.dist(result)
	attr(result, "method") = "dist"
	return(result) 
 }
######################################################################
# Function to perform the Breslow and Day (1980) test including
breslowday.test <- function(x) {
	#Find the common OR based on Mantel-Haenszel
	or.hat.mh <- mantelhaen.test(x)$estimate
	#Number of strata
	K <- dim(x)[3]
	#Value of the Statistic
	X2.HBD <- 0
	#Value of aj, tildeaj and Var.aj
	a <- tildea <- Var.a <- numeric(K)
	
	for (j in 1:K) {
		#Find marginals of table j
		mj <- apply(x[,,j], MARGIN=1, sum)
		nj <- apply(x[,,j], MARGIN=2, sum)
	
		#Solve for tilde(a)_j
		coef <- c(-mj[1]*nj[1] * or.hat.mh, nj[2]-mj[1]+or.hat.mh*(nj[1]+mj[1]),1-or.hat.mh)
		sols <- Re(polyroot(coef))
		#Take the root, which fulfills 0 < tilde(a)_j <= min(n1_j, m1_j)
		tildeaj <- sols[(0 < sols) &  (sols <= min(nj[1],mj[1]))]
		#Observed value
		aj <- x[1,1,j]
		
		#Determine other expected cell entries
		tildebj <- mj[1] - tildeaj
		tildecj <- nj[1] - tildeaj
		tildedj <- mj[2] - tildecj
	
		#Compute \hat{\Var}(a_j | \widehat{\OR}_MH)
		Var.aj <- (1/tildeaj + 1/tildebj + 1/tildecj + 1/tildedj)^(-1)
	
		#Compute contribution
		X2.HBD <- X2.HBD + as.numeric((aj - tildeaj)^2 / Var.aj)
	
		#Assign found value for later computations
		a[j] <- aj ;  tildea[j] <- tildeaj ; Var.a[j] <- Var.aj
	}
	
	#Compute Tarone corrected test
	X2.HBDT <-as.numeric( X2.HBD -  (sum(a) - sum(tildea))^2/sum(Var.aj) )
	
	#Compute p-value based on the Tarone corrected test
	p <- 1-pchisq(X2.HBDT, df=K-1)
	
	list(X2.HBD=X2.HBD,X2.HBDT=X2.HBDT,p=p)
}
