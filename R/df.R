df0=function(..., row.names = NULL, check.rows = FALSE,
           check.names = TRUE, fix.empty.names = TRUE,
           stringsAsFactors = FALSE){
	if(grepl('fix.empty.names',paste(format(args(data.frame)),collapse=''))){
	data.frame(..., row.names = row.names, check.rows = check.rows,
           check.names = check.names, fix.empty.names = fix.empty.names,
           stringsAsFactors = stringsAsFactors)
	}else{
	data.frame(..., row.names = row.names, check.rows = check.rows,
           check.names = check.names,
           stringsAsFactors = stringsAsFactors)
	}
}
df2mat=function(x,val,columnFactors,rowFactors,default=NA){
	x=as.data.frame(x)
	for(f in c(columnFactors,rowFactors)) x[,f]=factor(x[,f])
	if(! all(sapply(c(columnFactors,rowFactors),function(f) is.factor(x[,f])))) stop('All variables specified by columnFactors and rowFactors must be factors\n')
	cols=lapply(columnFactors,function(f) levels(x[,f]))
	names(cols)=NULL
	if(is.character(columnFactors)) names(cols)=columnFactors
	cLens=lapply(cols,length)
	nCol=do.call(prod,cLens)
	rows=lapply(rowFactors,function(f) levels(x[,f]))
	names(rows)=NULL
	if(is.character(rowFactors)) names(rows)=rowFactors
	rLens=lapply(rows,length)
	nRow=do.call(prod,rLens)
	if(length(columnFactors)>1){
		Cf=do.call(paste,x[,columnFactors])
		CnamesDf=lapply(1:(length(cols)-1),function(i) rep(cols[[i]],each=do.call(prod,cLens[-c(1:i)])))
		CnamesDf=c(CnamesDf,cols[length(cols)])
		CnamesDf=lapply(CnamesDf,rep,length.out=nCol)
		CnamesDf=as.data.frame(CnamesDf,stringsAsFactors=FALSE)
		colnames(CnamesDf)=names(cols)
		Cnames=do.call(paste,CnamesDf)
	}else{
		Cf=as.vector(x[,columnFactors])
		Cnames=levels(x[,columnFactors])
		CnamesDf=data.frame(Cnames,stringsAsFactors=FALSE)
		colnames(CnamesDf)=names(cols)
	}
	for(i in 1:length(columnFactors)) CnamesDf[,i]=factor(CnamesDf[,i],levels=cols[[i]])
	if(length(rowFactors)>1){
		Rf=do.call(paste,x[,rowFactors])
		RnamesDf=lapply(1:(length(rows)-1),function(i) rep(rows[[i]],each=do.call(prod,rLens[-c(1:i)])))
		RnamesDf=c(RnamesDf,rows[length(rows)])
		RnamesDf=lapply(RnamesDf,rep,length.out=nRow)
		RnamesDf=as.data.frame(RnamesDf,stringsAsFactors=FALSE)
		colnames(RnamesDf)=names(rows)
		Rnames=do.call(paste,RnamesDf)
	}else{
		Rf=as.vector(x[,rowFactors])
		Rnames=levels(x[,rowFactors])
		RnamesDf=data.frame(Rnames,stringsAsFactors=FALSE)
		colnames(RnamesDf)=names(rows)
	}
	for(i in 1:length(rowFactors)) RnamesDf[,i]=factor(RnamesDf[,i],levels=rows[[i]])
	mat=matrix(default,nrow=nRow,ncol=nCol)
	Cfi=match(Cf,Cnames)
	Rfi=match(Rf,Rnames)
	for(i in 1:nrow(x)) mat[Rfi[i],Cfi[i]]=x[i,val]
	attr(mat,'dimnames')=list(Rnames,Cnames)
	attr(mat,'factors')=list(rows=rows,cols=cols)
	attr(mat,'annotations')=list(RnamesDf=RnamesDf,CnamesDf=CnamesDf)
	mat
}
#convert list to two-column data.frame
list2df=function(x){
#output, a two-column df, with first column being the elements in x and the second column being the list entry name/id
	stopifnot(all(sapply(x,class) %in% c('numeric','integer','character')))
	n=sapply(x,length)
	data.frame(value=unlist(x),id=if(is.null(names(x))){rep(1:length(x),n)}else{rep(names(x),n)},stringsAsFactors=FALSE)
}
#rbind two data frame or matrices which may have different columns
rbind2m=function(x,y){
	if(is.null(x)) return(y)
	if(is.null(y)) return(x)
	xiny=colnames(x) %in% colnames(y)
	if(!all(xiny)){
		for(i in which(!xiny)){
			y=cbind(y,NA)
			colnames(y)[ncol(y)]=colnames(x)[i]
		}
	}
	yinx=colnames(y) %in% colnames(x)
	if(!all(yinx)){
		for(j in which(!yinx)){
			x=cbind(x,NA)
			colnames(x)[ncol(x)]=colnames(y)[j]
		}
	}
	rbind(x,y[,colnames(x)])
}
#cbind two data frame or matrices which may have different rows
cbind2m=function(x,y){
	if(is.null(x)) return(y)
	if(is.null(y)) return(x)
	xiny=rownames(x) %in% rownames(y)
	if(!all(xiny)){
		for(i in which(!xiny)){
			y=rbind(y,NA)
			rownames(y)[nrow(y)]=rownames(x)[i]
		}
	}
	yinx=rownames(y) %in% rownames(x)
	if(!all(yinx)){
		for(j in which(!yinx)){
			x=rbind(x,NA)
			rownames(x)[nrow(x)]=rownames(y)[j]
		}
	}
	cbind(x,y[rownames(x),])
}
