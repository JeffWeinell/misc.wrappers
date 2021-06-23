#' @title batch setup ustacks
#' 
#' This function creates a bash job file for each combination of input sequence read files (one per individual), values of m and M to run STACKS.
#' The output of this function is a character string vector with paths to bash job files, which can be piped to the function 'rsbatch' to submit the jobs.
#' 
#' @param ustacks.path Character string with path to the ustacks executable. Default NULL, in which case the path is derived from the default setting, which can be set with config_miscwrappers('<path/to/ustacks>').
#' @param files Input sequence data files, one for each individual.
#' @param save.in Character string with directory where output directories and files should be saved.
#' @param mMN Three-column data frame or matrix with m, M, and N values in the first, second, and third columns, respectively. Each row indicates a combination of m, M, and N to use for STACKS. Default cbind(m=rep(3,8),M=1:8,N=3:10).
#' m = Minimum depth of coverage required to create a stack. Default 3.
#' M = Maximum distance (in nucleotides) allowed between stacks. Default 1-8.
#' N = Maximum distance allowed to align secondary reads to primary stacks. Default 'M+2'.
#' @param N Maximum distance allowed to align secondary reads to primary stacks. Default 'M+2'.
#' @param p Number of threads to use for job parallization. Default 16.
#' @return Character string vector with the paths to bash job files.
#' @export ustacks_setup
ustacks_setup <- function(files,save.in,mMN=cbind(m=rep(3,8),M=1:8,N=3:10),p=16,ustacks.path=NULL){
	if(is.null(ustacks.path)){
		if(any(config_miscwrappers()[,"program"]=="ustacks")){
			ustacks.path <- config_miscwrappers()[config_miscwrappers()[,"program"]=="ustacks","exe_path"]
		} else {
			stop("Set 'ustacks.path' argument, or use config_miscwrappers(exe.path=<path/to/ustacks>) to set the package default path for ustacks")
		}
	}
	mMN            <- data.frame(m=mMN[,1],M=mMN[,2],N=mMN[,3])
	outdirs        <- file.path(save.in, paste0("m",mMN$m,"_M",mMN$M,"_N",mMN$N))
	make.outdirs   <- sapply(paste("mkdir -p",outdirs),system)
	all.df         <- data.frame(s=ustacks.path, f=rep(files,each=nrow(mMN)), i=rep(formatC(1:length(files),width=4,format="d",flag="0"),each=nrow(mMN)),o=rep(outdirs,length(files)), m=rep(mMN$m, length(files)), M=rep(mMN$M, length(files)), N=rep(mMN$N, length(files)), p=p)
	stacks.strings <- sapply(X=1:nrow(all.df),FUN=function(x) {do.call(sprintf,unname(c("%s -f %s -i %s -o %s -m %s -M %s -N %s -p %s",as.list(all.df[x,]))))})
	sh.lines       <- lapply(1:length(stacks.strings),FUN=function(x){c("#!/bin/bash",stacks.strings[x])})
	sh.paths       <- file.path(save.in,paste0("ustacks",1:length(stacks.strings),".sh"))
	sh.write       <- lapply(1:length(sh.lines),FUN=function(x){writeLines(text=sh.lines[[x]],con=sh.paths[x])})
	sh.paths
}

#' @title cleanup ustacks
#' 
#' Deletes ustacks bash scripts in a directory. Run this after all ustacks jobs are completed.
#' 
#' @param dir Character string with path to directory containing ustacks bash scripts
#' @return NULL
#' @export ustacks_cleanup
ustacks_cleanup <- function(dir){
	ushfiles <- list.files(dir,full.names=T,pattern="^ustacks[0-9]+.sh$")
	file.remove(ushfiles)
}

#' @title batch setup cstacks
#' 
#' For each combination of N and set of ustacks, this function creates a bash job file to run cstacks.
#' The output of this function is a character string vector with paths to bash job files, which can be piped to the function 'rsbatch' to submit the jobs.
#' 
#' @param cstacks.path Character string with path to the cstacks executable. Default NULL, in which case the path is derived from the default setting, which can be set with config_miscwrappers('<path/to/cstacks>').
#' @param save.in Character string with directory where output directories and files should be saved. This should be the same as the value used for the 'save.in' argument of ustacks_setup.
#' @param n Number or numerical vector with number of mismatches allowed between sample loci when build the catalog. Default 1.
#' @param p Number of threads to use for job parallization. Default 16.
#' @return Character string vector with the paths to bash job files.
#' @export cstacks_setup
cstacks_setup <- function(files,save.in,n=1,p=16,cstacks.path=NULL){
	if(is.null(cstacks.path)){
		if(any(config_miscwrappers()[,"program"]=="cstacks")){
			cstacks.path <- config_miscwrappers()[config_miscwrappers()[,"program"]=="cstacks","exe_path"]
		} else {
			stop("Set 'cstacks.path' argument, or use config_miscwrappers(exe.path=<path/to/cstacks>) to set the package default path for cstacks")
		}
	}
	#cstacks -o ./stacks -s ./stacks/f0_male -s ./stacks/f0_female -n 4 -p 15
	mMdirs         <- list.dirs(save.in,recursive=F)
	outdirs        <- paste0(rep(mMdirs,each=length(n)),rep(paste0("_n",n),length(mMdirs)))
	if(lenth(n)==1){
		file.rename(mMdirs,outdirs)
	} else {
		make.outdirs <- sapply(paste("mkdir -p",outdirs),system)
		outdirs.mat  <- matrix(data=outdirs,ncol=length(n),byrow=T)
		copy.ufiles  <- lapply(1:nrow(outdirs.mat),FUN=function(x){sapply(1:length(n),FUN=function(y){file.copy(list.files(mMdirs[x],full.names=T), outdirs.mat[x,y])})})
		delete.old   <- sapply(paste("rm -R",mMdirs),system)
	}
	samples        <- unlist(lapply(X=1:length(outdirs), FUN=function(x) {paste(paste("-s",sq(gsub(".snps.tsv.gz$","",list.files(outdirs[x],full.names=T,pattern=".snps.tsv.gz$")))),collapse=" ")}))
	all.df         <- data.frame(c=sq(cstacks.path),o=paste("-o",sq(outdirs)),s=samples,n=paste("-n",rep(n,length(mMdirs))),p=paste("-p",p))
	stacks.strings <- sapply(1:nrow(all.df),FUN=function(x) {paste(all.df[x,],collapse=" ")})
	sh.lines       <- lapply(1:length(stacks.strings),FUN=function(x){c("#!/bin/bash",stacks.strings[x])})
	sh.paths       <- file.path(save.in,paste0("cstacks",1:length(stacks.strings),".sh"))
	sh.write       <- lapply(1:length(sh.lines),FUN=function(x){writeLines(text=sh.lines[[x]],con=sh.paths[x])})
	sh.paths
}

#' @title cleanup cstacks
#' 
#' Deletes cstacks bash scripts in a directory. Run this after all cstacks jobs are completed.
#' 
#' @param dir Character string with path to directory containing cstacks bash scripts
#' @return NULL
#' @export ustacks_cleanup
ustacks_cleanup <- function(dir){
	cshfiles <- list.files(dir,full.names=T,pattern="^cstacks[0-9]+.sh$")
	file.remove(cshfiles)
}

save.in <- "/panfs/pfs.local/home/j926w878/scratch/scratch_v2/ddRAD/StacksFiles/Ahaetulla_AllPops"



