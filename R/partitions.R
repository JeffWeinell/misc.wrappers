#' @title vcf2parts
#' 
#' Generate a phylip format partition file from VCF file.
#' 
#' @param x 'vcfR' object (see package::vcfR) or character string with path to VCF file containing nucleotide data.
#' @param save.as Path where partition file should be saved.
#' @return Character string with path to partition file.
#' @export vcf2parts
vcf2parts <- function(x,save.as){
	if(file.exists(save.as)){
		stop("Output file already exists. Use a different name for 'save.as' argument.")
	}
	if(is(x,"vcfR")){
		vcf.obj <- vcf <- x
	} else {
		vcf <- x
		vcf.obj <- vcfR::read.vcfR(vcf,verbose=F,checkFile=F)
	}
	gt    <- vcf.obj@gt
	fx    <- vcf.obj@fix
	met   <- vcf.obj@meta
	samplenames <- colnames(gt)[-1]
	loci        <- c(fx[,"CHROM"])
	loci.unique <- unique(loci)
	# list of numerical vectors indicating range of sites for each locus (start and end)
	ranges.list     <- lapply(X=1:length(loci.unique),FUN=function(x){range(which(loci == loci.unique[x]))})
	# two-column matrix with range (start and end) of each locus
	ranges.mat      <- do.call(rbind,ranges.list)
	partition.lines <- paste0("DNA, part",1:nrow(ranges.mat)," = ",ranges.mat[,1],"-",ranges.mat[,2])
	writeLines(text=partition.lines,con=save.as)
	save.as
}

#' @title VCF format to Phylip format
#' 
#' Converts VCF to phylip and saves the result to file. Output can be sequential or interleaved.
#' Only use for diploids. The exported sequence is the haplotype consensus for each individual.
#' 
#' @param x 'vcfR' object (see package::vcfR) or character string with path to VCF file containing nucleotide data.
#' @param save.as Path where phylip alignment should be saved.
#' @param missing Character to use for missing data in the output file. Default "-".
#' @param sequential Whether or not data should be written on a single line (TRUE) or in multiple blocks (interleaced format) with the number of bases per line controlled by 'width'.
#' @param width NULL or a number specifying the number of bases to write per line when output is interleaved.
#' @param padding Number of spaces between the longest sample name and the first base.
#' @return The value of 'save.as'.
#' @export vcf2phy
vcf2phy <- function(x,save.as,missing="-",sequential=TRUE,width=300,padding=5,partition.file=TRUE){
	if(file.exists(save.as)){
		stop("Output file already exists. Use a different name for 'save.as' argument.")
	}
	if(is(x,"vcfR")){
		vcf.obj <- vcf <- x
	} else {
		vcf <- x
		vcf.obj <- vcfR::read.vcfR(vcf,verbose=F,checkFile=F)
	}
	gt    <- vcf.obj@gt
	fx    <- vcf.obj@fix
	met   <- vcf.obj@meta
	RA    <- cbind(fx[,"REF"],fx[,"ALT"])
	ra  <- paste0(fx[,"REF"],fx[,"ALT"])
	rat <- table(ra)
	gt0 <- gsub(":.+","",t(gt[,-1]))
	for(i in 1:length(rat)){
		ra.temp   <- names(rat)[i]
		RA.temp   <- unlist(strsplit(ra.temp,split=""))
		h.temp    <- toupper(seqinr::bma(RA.temp))
		RHA.temp  <- c(RA.temp[1],h.temp,RA.temp[2])
		cols.temp <- which(ra==ra.temp)
		sub.temp  <- gt0[,cols.temp]
		sub.temp1 <- gsub("0[/,|]0",RHA.temp[1],sub.temp)
		sub.temp2 <- gsub("0[/,|]1",RHA.temp[2],sub.temp1)
		sub.temp3 <- gsub("1[/,|]1",RHA.temp[3],sub.temp2)
		gt0[,cols.temp] <- sub.temp3
	}
	gt1 <- gsub("./.",missing,gt0,fixed=TRUE)
	gt1 <- gsub(".|.",missing,gt1,fixed=TRUE)
	gt1[is.na(gt1)] <- missing
	samplenames <- colnames(gt)[-1]
	nsites      <- nrow(gt)
	nindv       <- length(samplenames)
	npadcol <- (max(nchar(samplenames))-nchar(samplenames))+padding
	padcol  <- do.call(rbind,lapply(X=npadcol,FUN=function(x){paste0(rep(" ",x),collapse="")}))
	usenames <- paste0(samplenames,padcol)
	rownames(gt1) <- usenames
	if(sequential){
		gt2 <- do.call(rbind,as.list(apply(X=gt1,MARGIN=1,FUN=function(x){paste0(x,collapse="")})))
		gt3 <- cbind(rownames(gt1),gt2)
		gt4 <- unname(do.call(rbind,as.list(apply(X=gt3,MARGIN=1,FUN=function(x){paste0(x,collapse="")}))))
		gt5 <- c(paste(nindv,nsites),gt4)
		writeLines(text=gt5,con=save.as)
	} else {
		ranges.out <- matrix(data=c(1:nsites,rep(NA,((ceiling(nsites/width)*width)-nsites))),ncol=width,byrow=TRUE)
		ranges     <- do.call(rbind,lapply(X=1:nrow(ranges.out),FUN=function(x){range(ranges.out[x,],na.rm=TRUE)}))
		gt2 <- lapply(X=1:nrow(ranges),FUN=function(r){do.call(rbind,as.list(apply(X=gt1[,c(ranges[r,1]:ranges[r,2])], MARGIN=1,FUN=function(x){paste0(x,collapse="")})))})
		gt3 <- lapply(X=1:length(gt2),FUN=function(x){cbind(usenames,gt2[[x]])})
		gt4 <- lapply(X=gt3,FUN=function(j){c(unname(do.call(c,as.list(apply(X=j,MARGIN=1,FUN=function(x){paste0(x,collapse="")})))),"")})
		gt5 <- c(paste(nindv,nsites),do.call(c,gt4))
		writeLines(text=gt5,con=save.as)
	}
	if(partition.file){
		vcf2parts(x=x,save.as=paste0(save.as,".parts"))
	}
	save.as
}

