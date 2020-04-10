#wrapper to read/write tab delimited files
read.tsv=function(file, as.is=TRUE, header = TRUE, sep = "\t", quote = "\"",
           dec = ".", fill = TRUE, comment.char = "", verbose=TRUE, ...){
	if(verbose) cat('Read file ',file,'\n',sep='')
	read.delim(file,as.is=as.is,header=header,sep=sep,quote=quote,dec=dec,fill=fill,comment.char=comment.char,...)
}
write.tsv=function(x,file="",row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE,verbose=TRUE){
	write.table(x,file=file,row.names=row.names,col.names=col.names,sep=sep,quote=quote)
	if(verbose) cat('Written to file ',file,'\n',sep='')
}
save.rds=function(object, file = "", ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL){
	cat('saveRDS to file',file,'\n')
	saveRDS(object=object, file = file, ascii = ascii, version = version, compress = compress, refhook = refhook)
}
read.rds=function(file, refhook = NULL){
	cat('readRDS from file',file,'\n')
	readRDS(file=file, refhook = refhook)
}
writeCsv=function(x,file="",row.names=FALSE,quote=TRUE,verbose=TRUE){
	write.csv(x,file=file,row.names=row.names,quote=quote)
	if(verbose) cat('Written to file ',file,'\n',sep='')
}
fread2=function(file,verbose=TRUE,showProgress=FALSE,...){
	if(verbose) cat('Read file ',file,'\n',sep='')
	data.table::fread(file,verbose=!verbose,showProgress=showProgress,...)
}
makeDir=function(path, showWarnings = TRUE, recursive = TRUE, mode = "0777"){
	path=sub('\\/$','',path)
	if(! file.exists(path)) dir.create(path,recursive=recursive)
	return(invisible())
}
read.zip <- function(file, FUN=read.table, ...) {
  zipFileInfo <- unzip(file, list=TRUE)
  if(nrow(zipFileInfo) > 1)
    stop("More than one data file inside zip")
  else
    FUN(unz(file, as.character(zipFileInfo$Name)), ...)
}
#R functions to write and read matrix data in binary format
writeBinMatrix<-function(x,filename){
	nr=nrow(x)
	nc=ncol(x)
	conn=file(filename,'wb')
	ntmp = writeBin(c(nr,nc),conn, size = 8, endian = "little",useBytes=TRUE)
	ntmp = writeBin(rownames(x),conn, size = 8, endian = "little",useBytes=TRUE)
	ntmp = writeBin(colnames(x),conn, size = 8, endian = "little",useBytes=TRUE)
	for(i in 1:nc){
		ntmp = writeBin(x[,i],conn, size = 8, endian = "little",useBytes=TRUE)
	}
	close(conn)
	return(invisible())
}
readBinMatrix<-function(filename){ #to read matrix from a binary file
	conn=file(filename,'rb')
	nr=readBin(conn, integer(), size = 8, n = 1, endian = "little")
	nc=readBin(conn, integer(), size = 8, n = 1, endian = "little")
	dh=matrix(NA,nrow=nr,ncol=nc)
	rownames(dh)=readBin(conn, character(), size = 8, n = nr, endian = "little")
	colnames(dh)=readBin(conn, character(), size = 8, n = nc, endian = "little")
	for(i in 1:nc){
		dh[,i] = readBin(conn, numeric(), size = 8, n = nr, endian = "little")
	}
	close(conn)
	return(dh)
}
