#' @title batch setup ustacks
#' 
#' This function creates a bash job file for each combination of input sequence read files (one per individual), values of m and M to run STACKS.
#' The output of this function is a character string vector with paths to bash job files, which can be piped to the function 'rsbatch' to submit the jobs.
#' 
#' @param ustacks.path Character string with path to the ustacks executable.
#' @param files Input sequence data files, one for each individual.
#' @param save.in Character string with directory where output directories and files should be saved.
#' @param mM Two-column data frame or matrix with m and M values in the first and second columns, respectively. Each row indicates a combination of m and M to use for STACKS
#' @param p Number of threads to use for job parallization. Default 16.
#' @return Character string vector with the paths to bash job files.
#' @export ustacks_setup
ustacks_setup <- function(ustacks.path,files,save.in,mM,p=16){
	mM             <- data.frame(m=mM[,1],M=mM[,2])
	outdirs        <- file.path(save.in, paste0("m",mM$m,"_M",mM$M))
	make.outdirs   <- sapply(paste("mkdir -p",outdirs),system)
	all.df         <- data.frame(s=ustacks.path, f=rep(files,each=nrow(mM)), i=rep(formatC(1:length(files),width=4,format="d",flag="0"),nrow(mM)),o=rep(outdirs,length(files)), m=rep(mM$m, length(files)), M=rep(mM$M, length(files)), p=p)
	stacks.strings <- sapply(X=1:nrow(all.df),FUN=function(x) {do.call(sprintf,unname(c("%s -f %s -i %s -o %s -m %s -M %s -p %s",as.list(all.df[x,]))))})
	sh.lines       <- lapply(1:length(stacks.strings),FUN=function(x){c("#!/bin/bash",stacks.strings[x])})
	sh.paths       <- file.path(save.in,paste0("ustacks",1:length(stacks.strings),".sh"))
	sh.write       <- lapply(1:length(sh.lines),FUN=function(x){writeLines(text=sh.lines[[x]],con=sh.paths[x])})
	sh.paths
}
