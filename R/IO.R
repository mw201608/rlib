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
	if(! file.exists(path)) dir.create(path,recursive=recursive, mode = mode)
	return(invisible())
}
read.zip <- function(file, FUN=read.table, ...) {
  zipFileInfo <- unzip(file, list=TRUE)
  if(nrow(zipFileInfo) > 1)
    stop("More than one data file inside zip\n")
  else
    FUN(unz(file, as.character(zipFileInfo$Name)), ...)
}
read.gzip <- function(file, FUN = read.delim, ..., verbose = TRUE) {
    if(verbose) cat('Reading from', file, "\n")
	if(! grepl('^(ftp|http)', file, ignore.case = TRUE)){
		con <- gzfile(file)
	}else{
		con <- gzcon(url(file))
	}
	on.exit(close(con))
	if(deparse(substitute(FUN)) == "scan") return(scan(con, ...))
	txt <- readLines(con)
 	return(FUN(textConnection(txt), ...))
}
#R functions to write and read numeric matrix data in binary format
writeBinMatrix <- function(x, filename, dimBytes = 8, charBytes = 20, dataBytes = 8, endian = "little", useBytes = TRUE){
	nr = nrow(x)
	nc = ncol(x)
	if(charBytes > 0){
		if(is.null(rownames(x)) || is.null(colnames(x))) stop('Data must have both rownames and columns when charBytes is positive.\n')
	}
	conn = file(filename, 'wb')
	ntmp = writeBin(c(nr,nc), conn, size = dimBytes, endian = endian, useBytes = useBytes)
	if(charBytes > 0){
		ntmp = writeBin(colnames(x), conn, size = charBytes, endian = endian, useBytes = useBytes)
		ntmp = writeBin(rownames(x), conn, size = charBytes, endian = endian, useBytes = useBytes)
	}
	for(i in 1:nc){
		ntmp = writeBin(x[, i], conn, size = dataBytes, endian = endian, useBytes = useBytes)
	}
	close(conn)
	return(invisible())
}
readBinMatrix <- function(filename, dimBytes = 8, charBytes = 20, dataBytes = 8, endian = "little"){ #to read matrix from a binary file created by writeBinMatrix
	conn = file(filename, 'rb')
	on.exit(close(conn), add = FALSE, after = TRUE)
	nr = readBin(conn, integer(), size = dimBytes, n = 1, endian = endian)
	nc = readBin(conn, integer(), size = dimBytes, n = 1, endian = endian)
	if(charBytes > 0){
		cn = readBin(conn, character(), size = charBytes, n = nc, endian = endian)
		rn = readBin(conn, character(), size = charBytes, n = nr, endian = endian)
	}
	if(dataBytes > 0){
		dh = matrix(NA, nrow=nr, ncol=nc)
		for(i in 1:nc){
			dh[, i] = readBin(conn, numeric(), size = dataBytes, n = nr, endian = endian)
		}
		rownames(dh) <- rn
		colnames(dh) <- cn
	}else{
		if(charBytes > 0) return(list(rownames = rn, colnames = cn))
		return(c(nr, nc))
	}
	return(dh)
}
readSingleColumnFromBinMatrix <- function (j, id, filename, dimBytes = 8, charBytes = 20, dataBytes = 8, endian = "little", verbose = FALSE){
    if(verbose) cat('Read', filename, '...\n')
	if(! file.exists(filename)) return(NULL)
	if(missing(j)){
		if(missing(id)) return(readBinMatrix(filename = filename, dimBytes = dimBytes, charBytes = charBytes, dataBytes = dataBytes, endian = endian))
		if(charBytes <= 0) stop("charBytes must be a valid positive number in order to read column names\n")
	}else{
		stopifnot(j > 0)
	}
	conn = file(filename, "rb")
	on.exit(close(conn), add = FALSE, after = TRUE)
    nr = readBin(conn, integer(), size = dimBytes, n = 1, endian = endian)
    nc = readBin(conn, integer(), size = dimBytes, n = 1, endian = endian)
	if((! missing(j)) && j > nc) stop('Value j is not valid\n')
    if(charBytes > 0) colnam = readBin(conn, character(), size = charBytes, n = nc, endian = endian)
	if(missing(j)){
		j <- which(colnam == id)
		if(length(j) == 0) return(NULL)
	}
    if(charBytes > 0) rownam = readBin(conn, character(), size = charBytes, n = nr, endian = endian)
	if(j > 1) np = seek(conn, where = dataBytes * (j-1) * nr, origin = "current", rw = "r")
    dh = readBin(conn, numeric(), size = dataBytes, n = nr, endian = endian)
    if(charBytes > 0) dh = setNames(dh, rownam)
	return(dh)
}
computeHeaderBytes = function(Dimnames, dataBytes=8){
	2*dataBytes + sum(sapply(Dimnames[[1]], function(x) nchar(x) + 1)) + sum(sapply(Dimnames[[2]], function(x) nchar(x) + 1))
}
readDimsFromRemoteBinMatrix = function(url, dimBytes = 8, charBytes = 20, endian = 'little', dims = NULL, verbose = FALSE){
	if(is.null(dims)){
		r = paste0("bytes=0-", (2*dimBytes - 1))
		if(verbose) cat('Range', r, '\n')
		res = GET(url=url, config = add_headers(Range=r))
		if(http_status(res)$category != "Success") print_err_and_return(ifelse(verbose, 'Dimension retrieval failed\n', NULL))(NULL)
		dims = readBin(res$content, what='integer', size=dimBytes, n=2, endian = endian)
		if(charBytes == 0) return(dims)
	}
	r = paste0("bytes=", 2 * dimBytes, "-", (2 * dimBytes + (dims[2] + dims[1]) * charBytes - 1))
	if(verbose) cat('Range', r, '\n')
	res = GET(url=url, config = add_headers(Range=r))
	if(http_status(res)$category != "Success") print_err_and_return(ifelse(verbose, 'Dimnames retrieval failed\n', NULL))(NULL)
	nam = readBin(res$content, character(), size = charBytes, n = dims[1] + dims[2], endian = endian)
	return(list(rownames = nam[(1 + dims[2]) : (dims[1] + dims[2])], colnames = nam[1:dims[2]]))
}
readSingleColumnFromRemoteBinMatrix = function(j, id, url, dimBytes = 8, charBytes = 20, dataBytes = 8, endian = 'little', dims = NULL, verbose = FALSE){
	if(verbose) cat('Read', url, '...\n')
	Dimnames <- readDimsFromRemoteBinMatrix(url = url, dimBytes = dimBytes, charBytes = charBytes, endian = endian, dims = dims, verbose = verbose)
	if(is.null(Dimnames)) return(NULL)
	if(is.null(dims)) dims <- c(length(Dimnames$rownames), length(Dimnames$colnames))
	if(missing(j)){
		if(missing(id)) print_err_and_return(ifelse(verbose, 'Either j or id must be provided\n', NULL))(NULL)
		j <- which(Dimnames$colnames == id)
	}
	if(length(j) == 0) return(NULL)
	if(j < 1 || j > dims[2]) return(NULL)
	#
	totalHeaderBytes = computeHeaderBytes(Dimnames = Dimnames, dataBytes = dimBytes)
	b1 = totalHeaderBytes + (j-1) * dims[1] * dataBytes
	b2 = b1 + dims[1]*dataBytes
	r = paste0("bytes=", b1, "-", b2-1)
	if(verbose) cat('Range', r, '\n')
	res = GET(url=url, config=add_headers(Range=r))
	if(http_status(res)$category != "Success") print_err_and_return(ifelse(verbose, 'Data retrieval failed\n', NULL))(NULL)
	readBin(res$content, numeric(), size = dataBytes, n = dims[1], endian = endian)
}
