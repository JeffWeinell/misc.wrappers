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
	mMN   <- data.frame(m=mMN[,1],M=mMN[,2],N=mMN[,3])
	outdirs  <- file.path(save.in, paste0("m",mMN$m,"_M",mMN$M,"_N",mMN$N))
	#outdirs  <- file.path(save.in, paste0("m",mMN$m,"_M",mMN$M))
	make.outdirs <- sapply(paste("mkdir -p",outdirs),system)
	all.df   <- data.frame(s=ustacks.path, f=rep(files,each=nrow(mMN)), i=rep(formatC(1:length(files),width=4,format="d",flag="0"),each=nrow(mMN)),o=rep(outdirs,length(files)), m=rep(mMN$m, length(files)), M=rep(mMN$M, length(files)), N=rep(mMN$N, length(files)), p=p)
	stacks.strings <- sapply(X=1:nrow(all.df),FUN=function(x) {do.call(sprintf,unname(c("%s -f %s -i %s -o %s -m %s -M %s -N %s -p %s",as.list(all.df[x,]))))})
	sh.lines  <- lapply(1:length(stacks.strings),FUN=function(x){c("#!/bin/bash",stacks.strings[x])})
	sh.paths  <- file.path(save.in,paste0("ustacks",1:length(stacks.strings),".sh"))
	sh.write  <- lapply(1:length(sh.lines),FUN=function(x){writeLines(text=sh.lines[[x]],con=sh.paths[x])})
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
cstacks_setup <- function(save.in,n=1,p=16,cstacks.path=NULL){
	if(is.null(cstacks.path)){
		if(any(config_miscwrappers()[,"program"]=="cstacks")){
			cstacks.path <- config_miscwrappers()[config_miscwrappers()[,"program"]=="cstacks","exe_path"]
		} else {
			stop("Set 'cstacks.path' argument, or use config_miscwrappers(exe.path=<path/to/cstacks>) to set the package default path for cstacks")
		}
	}
	#mMNdirs  <- list.dirs(save.in,recursive=F)
	mMNdirs  <- list.files(save.in,recursive=F,pattern="^m[0-9]+_M[0-9]+_N[0-9]+$",include.dirs=T,full.names=T)
	if(length(mMNdirs)==0){
		mMNndirs <- list.files(save.in,recursive=F,pattern="^m[0-9]+_M[0-9]+_N[0-9]+_n[0-9]+$",include.dirs=T,full.names=T)
		if(length(mMNndirs)>0){
			setup.dirs <- FALSE
			outdirs <- mMNndirs
			mMNdirs <- gsub("_n[0-9]+$","",mMNndirs)
		} else {
			stop("no input directories found")
		}
	} else {
		outdirs <- paste0(rep(mMNdirs,each=length(n)),rep(paste0("_n",n),length(mMNdirs)))
		if(length(n)==1){
			file.rename(mMNdirs,outdirs)
		} else {
			make.outdirs <- sapply(paste("mkdir -p",outdirs),system)
			outdirs.mat  <- matrix(data=outdirs,ncol=length(n),byrow=T)
			copy.ufiles  <- lapply(1:nrow(outdirs.mat),FUN=function(x){sapply(1:length(n),FUN=function(y){file.copy(list.files(mMNdirs[x],full.names=T), outdirs.mat[x,y])})})
			#delete.old   <- sapply(paste("rm -R",mMNdirs),system)
		}
	}
	samples  <- unlist(lapply(X=1:length(outdirs), FUN=function(x) {paste(paste("-s",sq(gsub(".snps.tsv.gz$","",list.files(outdirs[x],full.names=T,pattern=".snps.tsv.gz$")))),collapse=" ")}))
	all.df   <- data.frame(c=sq(cstacks.path),o=paste("-o",sq(outdirs)),s=samples,n=paste("-n",rep(n,length(mMNdirs))),p=paste("-p",p))
	stacks.strings <- sapply(1:nrow(all.df),FUN=function(x) {paste(all.df[x,],collapse=" ")})
	sh.lines  <- lapply(1:length(stacks.strings),FUN=function(x){c("#!/bin/bash",stacks.strings[x])})
	sh.paths  <- file.path(save.in,paste0("cstacks",1:length(stacks.strings),".sh"))
	sh.write  <- lapply(1:length(sh.lines),FUN=function(x){writeLines(text=sh.lines[[x]],con=sh.paths[x])})
	sh.paths
}

#' @title cleanup cstacks
#' 
#' Deletes cstacks bash scripts in a directory. Run this after all cstacks jobs are completed.
#' 
#' @param dir Character string with path to directory containing cstacks bash scripts
#' @return NULL
#' @export cstacks_cleanup
cstacks_cleanup <- function(dir){
	cshfiles <- list.files(dir,full.names=T,pattern="^cstacks[0-9]+.sh$")
	file.remove(cshfiles)
}


#' @title batch setup sstacks
#' 
#' For each set samples (processed with ustacks and cstacks) and a population map, this function creates a bash job file to run sstacks.
#' The output of this function is a character string vector with paths to bash job files, which can be piped to the function 'rsbatch' to submit the jobs.
#' 
#' @param save.in Character string with directory where output directories and files should be saved. This should be the same as the value used for the 'save.in' argument of ustacks_setup.
#' @param popmap Character string to population map. Default NULL, in which case all individuals in 'save.in' are assigned to the same 'population'.
#' @param p Number of threads to use for job parallization. Default 16.
#' @param sstacks.path Character string with path to the sstacks executable. Default NULL, in which case the path is derived from the default setting, which can be set with config_miscwrappers('<path/to/sstacks>').
#' @return Character string vector with the paths to bash job files.
#' @export sstacks_setup
sstacks_setup <- function(save.in,popmap=NULL,p=16,sstacks.path=NULL){
	if(is.null(sstacks.path)){
		if(any(config_miscwrappers()[,"program"]=="sstacks")){
			sstacks.path <- config_miscwrappers()[config_miscwrappers()[,"program"]=="sstacks","exe_path"]
		} else {
			stop("Set 'sstacks.path' argument, or use config_miscwrappers(exe.path=<path/to/sstacks>) to set the package default path for sstacks")
		}
	}
	mMNndirs <- list.dirs(save.in,recursive=F)
	outdirs <- mMNndirs
	if(is.null(popmap)){
		popmap.df <- data.frame(indv=basename(unique(gsub(".snps.tsv.gz","",list.files(mMNndirs[1],full.names=T,pattern=".snps.tsv.gz")))),pop=1)
		popmap.df2 <- popmap.df[popmap.df[,"indv"]!="catalog",]
		popmap.path <- file.path(save.in,"popmap.txt")
		write.table(popmap.df2,file=popmap.path,col.names=F,row.names=F,sep="\t",quote=F)
	} else {
		popmap.path <- popmap
	}
	#samples  <- unlist(lapply(X=1:length(outdirs), FUN=function(x) {paste(paste("-s",sq(gsub(".snps.tsv.gz$","",list.files(outdirs[x],full.names=T,pattern=".snps.tsv.gz$")))),collapse=" ")}))
	# sstacks -P ./stacks -M ./popmap -p 8
	all.df   <- data.frame(ss=sq(sstacks.path),P=paste("-P",sq(outdirs)),M=paste("-M",popmap.path),p=paste("-p",p))
	stacks.strings <- sapply(1:nrow(all.df),FUN=function(x) {paste(all.df[x,],collapse=" ")})
	sh.lines  <- lapply(1:length(stacks.strings),FUN=function(x){c("#!/bin/bash",stacks.strings[x])})
	sh.paths  <- file.path(save.in,paste0("sstacks",1:length(stacks.strings),".sh"))
	sh.write  <- lapply(1:length(sh.lines),FUN=function(x){writeLines(text=sh.lines[[x]],con=sh.paths[x])})
	sh.paths
}

#' @title cleanup sstacks
#' 
#' Deletes sstacks bash scripts in a directory. Run this after all sstacks jobs are completed.
#' 
#' @param dir Character string with path to directory containing sstacks bash scripts
#' @return NULL
#' @export sstacks_cleanup
sstacks_cleanup <- function(dir){
	sshfiles <- list.files(dir,full.names=T,pattern="^sstacks[0-9]+.sh$")
	file.remove(sshfiles)
}

#' @title batch setup tsv2bam
#' 
#' For each set samples (processed with ustacks, cstacks, and sstacks) and a population map, this function creates a bash job file to run tsv2bam.
#' The output of this function is a character string vector with paths to bash job files, which can be piped to the function 'rsbatch' to submit the jobs.
#' 
#' @param save.in Character string with directory where output directories and files should be saved. This should be the same as the value used for the 'save.in' argument of ustacks_setup.
#' @param popmap Character string to population map. Default NULL, in which case all individuals in 'save.in' are assigned to the same 'population'.
#' @param t Number of threads to use for job parallization. Default 16.
#' @param samples Character string vector with paths to input '.fq.gz' data files (i.e., those supplied to ustacks). Alternative to using 'samples.dir' argument. Default NULL. If 'samples' is also NULL, then data are assumed to be single rather than paired end reads.
#' @param samples.dir Character string with directory containing '.fq.gz' data files. Alternative to using 'samples' argument. Default NULL. If 'samples' is also NULL, then data are assumed to be single rather than paired end reads.
#' @param tsv2bam.path Character string with path to the tsv2bam executable. Default NULL, in which case the path is derived from the default setting, which can be set with config_miscwrappers('<path/to/tsv2bam>').
#' @return Character string vector with the paths to bash job files.
#' @export tsv2bam_setup
tsv2bam_setup <- function(save.in,popmap=NULL,samples=NULL,samples.dir=NULL,t=16,tsv2bam.path=NULL){
	if(is.null(tsv2bam.path)){
		if(any(config_miscwrappers()[,"program"]=="tsv2bam")){
			tsv2bam.path <- config_miscwrappers()[config_miscwrappers()[,"program"]=="tsv2bam","exe_path"]
		} else {
			stop("Set 'tsv2bam.path' argument, or use config_miscwrappers(exe.path=<path/to/tsv2bam>) to set the package default path for tsv2bam")
		}
	}
	mMNndirs <- list.dirs(save.in,recursive=F)
	outdirs <- mMNndirs
	if(is.null(popmap)){
		if("popmap.txt" %in% list.files(save.in)){
			popmap.path <- file.path(save.in,"popmap.txt")
		} else {
			popmap.df <- data.frame(indv=basename(unique(gsub(".snps.tsv.gz","",list.files(mMNndirs[1],full.names=T,pattern=".snps.tsv.gz")))),pop=1)
			popmap.df2 <- popmap.df[popmap.df[,"indv"]!="catalog",]
			popmap.path <- file.path(save.in,"popmap.txt")
			write.table(popmap.df2,file=popmap.path,col.names=F,row.names=F,sep="\t",quote=F)
		}
	} else {
		popmap.path <- popmap
	}
	# tsv2bam -P ./stacks/ -M ./popmap -R ./samplesdir -s sample1.fq.gz -s sample2.fq.gz-t 8
	if(!is.null(samples)){
		sample.paths <- paste("-s",samples)
	} else {
		sample.paths <- NA
	}
	if(!is.null(samples.dir)){
		samples.dir <- paste("-R",samples.dir)
	} else {
		samples.dir <- NA
	}
	all.df0   <- data.frame(t2b=sq(tsv2bam.path),P=paste("-P",sq(outdirs)),M=paste("-M",popmap.path),R=samples.dir,s=sample.paths,t=paste("-t",t))
	all.df <- all.df0[,!sapply(1:ncol(all.df0),FUN=function(x){all(is.na(all.df0[,x]))})]
	stacks.strings <- sapply(1:nrow(all.df),FUN=function(x) {paste(all.df[x,],collapse=" ")})
	sh.lines  <- lapply(1:length(stacks.strings),FUN=function(x){c("#!/bin/bash",stacks.strings[x])})
	sh.paths  <- file.path(save.in,paste0("tsv2bam",1:length(stacks.strings),".sh"))
	sh.write  <- lapply(1:length(sh.lines),FUN=function(x){writeLines(text=sh.lines[[x]],con=sh.paths[x])})
	sh.paths
}

#' @title cleanup tsv2bam
#' 
#' Deletes tsv2bam bash scripts in a directory. Run this after all tsv2bam jobs are completed.
#' 
#' @param dir Character string with path to directory containing tsv2bam bash scripts
#' @return NULL
#' @export tsv2bam_cleanup
tsv2bam_cleanup <- function(dir){
	t2bshfiles <- list.files(dir,full.names=T,pattern="^tsv2bam[0-9]+.sh$")
	file.remove(t2bshfiles)
}

#' @title batch setup gstacks
#' 
#' For each set samples (processed with ustacks, cstacks, sstacks, and tsv2bam) and a population map, this function creates a bash job file to run gstacks.
#' The output of this function is a character string vector with paths to bash job files, which can be piped to the function 'rsbatch' to submit the jobs.
#' 
#' @param save.in Character string with directory where output directories and files should be saved. This should be the same as the value used for the 'save.in' argument of ustacks_setup.
#' @param popmap Character string to population map. Default NULL, in which case all individuals in 'save.in' are assigned to the same 'population'.
#' @param t Number of threads to use for job parallization. Default 16.
#' @param gstacks.path Character string with path to the gstacks executable. Default NULL, in which case the path is derived from the default setting, which can be set with config_miscwrappers('<path/to/gstacks>').
#' @return Character string vector with the paths to bash job files.
#' @export gstacks_setup
gstacks_setup <- function(save.in,popmap=NULL,t=16,gstacks.path=NULL){
	if(is.null(gstacks.path)){
		if(any(config_miscwrappers()[,"program"]=="gstacks")){
			gstacks.path <- config_miscwrappers()[config_miscwrappers()[,"program"]=="gstacks","exe_path"]
		} else {
			stop("Set 'gstacks.path' argument, or use config_miscwrappers(exe.path=<path/to/gstacks>) to set the package default path for gstacks")
		}
	}
	mMNndirs <- list.dirs(save.in,recursive=F)
	outdirs <- mMNndirs
	if(is.null(popmap)){
		if("popmap.txt" %in% list.files(save.in)){
			popmap.path <- file.path(save.in,"popmap.txt")
		} else {
			popmap.df <- data.frame(indv=basename(unique(gsub(".snps.tsv.gz","",list.files(mMNndirs[1],full.names=T,pattern=".snps.tsv.gz")))),pop=1)
			popmap.df2 <- popmap.df[popmap.df[,"indv"]!="catalog",]
			popmap.path <- file.path(save.in,"popmap.txt")
			write.table(popmap.df2,file=popmap.path,col.names=F,row.names=F,sep="\t",quote=F)
		}
	} else {
		popmap.path <- popmap
	}
	all.df   <- data.frame(gs=sq(gstacks.path),P=paste("-P",sq(outdirs)),M=paste("-M",popmap.path),t=paste("-t",t))
	#all.df <- all.df0[,!sapply(1:ncol(all.df0),FUN=function(x){all(is.na(all.df0[,x]))})]
	stacks.strings <- sapply(1:nrow(all.df),FUN=function(x) {paste(all.df[x,],collapse=" ")})
	sh.lines  <- lapply(1:length(stacks.strings),FUN=function(x){c("#!/bin/bash",stacks.strings[x])})
	sh.paths  <- file.path(save.in,paste0("gstacks",1:length(stacks.strings),".sh"))
	sh.write  <- lapply(1:length(sh.lines),FUN=function(x){writeLines(text=sh.lines[[x]],con=sh.paths[x])})
	sh.paths
}

#' @title cleanup gstacks
#' 
#' Deletes gstacks bash scripts in a directory. Run this after all gstacks jobs are completed.
#' 
#' @param dir Character string with path to directory containing gstacks bash scripts
#' @return NULL
#' @export gstacks_cleanup
gstacks_cleanup <- function(dir){
	gshfiles <- list.files(dir,full.names=T,pattern="^gstacks[0-9]+.sh$")
	file.remove(gshfiles)
}

#' @title batch setup populations
#' 
#' For each set samples (processed with ustacks, cstacks, sstacks, and tsv2bam) and a population map, this function creates a bash job file to run populations.
#' The output of this function is a character string vector with paths to bash job files, which can be piped to the function 'rsbatch' to submit the jobs.
#' 
#' @param save.in Character string with directory where output directories and files should be saved. This should be the same as the value used for the 'save.in' argument of ustacks_setup.
#' @param use.popmap Logical indicating whether or not to assign individuals to different populations according to a population map. Default FALSE.
#' @param popmap Character string to population map. Default NULL, in which case all individuals in 'save.in' are assigned to the same 'population'.
#' @param include.out Character string vector indicating which output files should be generated. Default is to include all the files possible. See STACKS for explanation.
#' fasta-loci — output locus consensus sequences in FASTA format..
#' fasta-samples — output the sequences of the two haplotypes of each (diploid) sample, for each locus, in FASTA format.
#' vcf — output SNPs and haplotypes in Variant Call Format (VCF).
#' genepop — output results in GenePop format.
#' structure — output results in Structure format.
#' radpainter — output results in fineRADstructure/RADpainter format.
#' plink — output genotypes in PLINK format.
#' hzar — output genotypes in Hybrid Zone Analysis using R (HZAR) format.
#' phylip — output nucleotides that are fixed-within, and variant among populations in Phylip format for phylogenetic tree construction.
#' phylip-var — include variable sites in the phylip output encoded using IUPAC notation.
#' phylip-var-all — include all sequence as well as variable sites in the phylip output encoded using IUPAC notation.
#' treemix — output SNPs in a format useable for the TreeMix program (Pickrell and Pritchard).
#' no-hap-exports — omit haplotype outputs.
#' fasta-samples-raw — output all haplotypes observed in each sample, for each locus, in FASTA format.
#' @param t Number of threads to use for job parallization. Default 16.
#' @param populations.path Character string with path to the populations executable. Default NULL, in which case the path is derived from the default setting, which can be set with config_miscwrappers('<path/to/populations>').
#' @return Character string vector with the paths to bash job files.
#' @export populations_setup
populations_setup <- function(save.in,use.popmap=FALSE,popmap=NULL,include.out=c("fasta-loci","fasta-samples","vcf","genepop","structure","radpainter","plink","hzar","phylip","phylip-var","phylip-var-all","treemix","fasta-samples-raw"),t=16,populations.path=NULL){
	if(is.null(populations.path)){
		if(any(config_miscwrappers()[,"program"]=="populations")){
			populations.path <- config_miscwrappers()[config_miscwrappers()[,"program"]=="populations","exe_path"]
		} else {
			stop("Set 'populations.path' argument, or use config_miscwrappers(exe.path=<path/to/populations>) to set the package default path for populations")
		}
	}
	mMNndirs <- list.dirs(save.in,recursive=F)
	outdirs <- mMNndirs
	if(use.popmap){
		if(is.null(popmap)){
			if("popmap.txt" %in% list.files(save.in)){
				popmap.path <- file.path(save.in,"popmap.txt")
			} else {
				popmap.df <- data.frame(indv=basename(unique(gsub(".snps.tsv.gz","",list.files(mMNndirs[1],full.names=T,pattern=".snps.tsv.gz")))),pop=1)
				popmap.df2 <- popmap.df[popmap.df[,"indv"]!="catalog",]
				popmap.path <- file.path(save.in,"popmap.txt")
				write.table(popmap.df2,file=popmap.path,col.names=F,row.names=F,sep="\t",quote=F)
			}
		} else {
			popmap.path <- popmap
		}
		popmap.path <- paste("-M",popmap.path)
	} else {
		popmap.path <- NA
	}
	outformats  <- paste(paste0("--",include.out),collapse=" ")
	all.df0   <- data.frame(ps=sq(populations.path),P=paste("-P",sq(outdirs)),M=popmap.path,of=outformats,t=paste("-t",t))
	all.df <- all.df0[,!sapply(1:ncol(all.df0),FUN=function(x){all(is.na(all.df0[,x]))})]
	stacks.strings <- sapply(1:nrow(all.df),FUN=function(x) {paste(all.df[x,],collapse=" ")})
	sh.lines  <- lapply(1:length(stacks.strings),FUN=function(x){c("#!/bin/bash",stacks.strings[x])})
	sh.paths  <- file.path(save.in,paste0("populations",1:length(stacks.strings),".sh"))
	sh.write  <- lapply(1:length(sh.lines),FUN=function(x){writeLines(text=sh.lines[[x]],con=sh.paths[x])})
	sh.paths
}

#' @title cleanup populations
#' 
#' Deletes populations bash scripts in a directory. Run this after all populations jobs are completed.
#' 
#' @param dir Character string with path to directory containing populations bash scripts
#' @return NULL
#' @export populations_cleanup
populations_cleanup <- function(dir){
	pshfiles <- list.files(dir,full.names=T,pattern="^populations[0-9]+.sh$")
	file.remove(pshfiles)
}

#' @title cleanup stacks
#' 
#' Deletes cstacks, sstacks, tsv2bam, gstacks, and populations bash scripts in a directory. Run this after completing the entire STACKS pipeline.
#' 
#' @param dir Character string with path to directory containing populations bash scripts
#' @return NULL
#' @export stacks_cleanup
stacks_cleanup <- function(dir){
	cstacks_cleanup(dir)
	sstacks_cleanup(dir)
	tsv2bam_cleanup(dir)
	gstacks_cleanup(dir)
	populations_cleanup(dir)
}

#' @title copy STACKS SNPs VCF files using informative names
#' 
#' Copies each of the 'populations.snps.vcf' files generated by the 'populations' program of the STACKS pipeline, and uses a more informative filename. This function will only work if the 'populations_setup' function was used for running populations.
#' 
#' @param indir Character string with the parent directory where STACKS output files and directories are held for a particular group.
#' @param outdir Character string with path to directory where 'populations.snps.vcf' files should be copied. These will be renamed with the groupname stacks parameters.
#' @param groupname Character string with name of the group (usually a taxon). Default is NULL, in which case the groupname is the directory name supplied by indir,i.e., the result of basename(indir).
#' @return NULL
#' @export cpSNPsVCF
cpSNPsVCF <- function(indir,outdir,groupname=NULL){
	if(is.null(groupname)){
		groupname=basename(indir)
	}
	if(!dir.exists(outdir)){
		system(paste("mkdir -p",outdir))
	}
	infiles  <- list.files(indir, pattern="^populations.snps.vcf$", full.name=T, recursive=T)
	outfiles <- file.path(outdir,paste0(groupname,"_",basename(dirname(infiles)),"_AllSNPs.vcf"))
	commands <- paste0("cp '",infiles,"' '",outfiles,"'")
	sapply(commands,system)
	return(0)
}




#' @title summarize_stacks
#' 
#' For each set of samples processed with populations (STACKS) and with vcf file generated for snps, this function creates a bash job file to calculate some basic popgen statistics.
#' 
#' @param dir Directory containing multiple subdirectories, each of which contains the output files generated by populations (STACKS pipeline).
#' @param save Logical indicating if the output should be saved. Default TRUE.
#' @param use.popmap Logical indicating whether or not to assign individuals to different populations according to a population map. Default FALSE.
#' @param popmap Character string to population map. Default NULL, in which case all individuals in 'save.in' are assigned to the same 'population'.
#' @return data frame with basic stats
#' @export summarize_stacks
summarize_stacks <- function(dir, save=TRUE, use.popmap=FALSE, popmap=NULL){
	save.in <- dir
	mMNndirs <- list.dirs(save.in, recursive=F)
	if(use.popmap){
		if(is.null(popmap)){
			if("popmap.txt" %in% list.files(save.in)){
				popmap.path <- file.path(save.in,"popmap.txt")
			} else {
				popmap.df <- data.frame(indv=basename(unique(gsub(".snps.tsv.gz","",list.files(mMNndirs[1],full.names=T,pattern=".snps.tsv.gz")))),pop=1)
				popmap.df2 <- popmap.df[popmap.df[,"indv"]!="catalog",]
				popmap.path <- file.path(save.in,"popmap.txt")
				write.table(popmap.df2,file=popmap.path,col.names=F,row.names=F,sep="\t",quote=F)
			}
		} else {
			popmap.path <- popmap
		}
		popmap.path <- paste("-M",popmap.path)
	} else {
		popmap.path <- NA
	}
	res <- data.frame()
	for(i in 1:length(mMNndirs)){
		vcf <- file.path(mMNndirs[i],"populations.snps.vcf")
		if(!file.exists(vcf)){
			stop("could not find 'populations.snps.vcf'")
		}
		vcf.obj    <- vcfR::read.vcfR(vcf,verbose=F,checkFile=F)
		gt         <- gsub(":.+","",vcf.obj@gt[,-1])
		gt.counts  <- table(gt,useNA='always')
		#counts.sites <- sapply(X=1:nrow(gt), FUN=function(x) {A=table(gt[x,], useNA='always'); B=as.numeric(A[is.na(names(A))])/sum(A); B})
		counts.sites <- lapply(X=1:nrow(gt), FUN=function(x) {table(gt[x,], useNA='always')})
		#counts.indv <- sapply(X=1:ncol(gt), FUN=function(x) {A=table(gt[,x], useNA='always'); B=as.numeric(A[is.na(names(A))])/sum(A); B})
		counts.indv <- lapply(X=1:ncol(gt), FUN=function(x) {table(gt[,x], useNA='always')})
		missingness.sites <- round(sapply(X=1:length(counts.sites),FUN=function(x) {A=counts.sites[[x]]; as.numeric(A[names(A) %in% c("./.",NA)])/sum(A)}),digits=4)
		missingness.indv <- round(sapply(X=1:length(counts.indv),FUN=function(x) {A=counts.indv[[x]]; as.numeric(A[names(A) %in% c("./.",NA)])/sum(A)}),digits=4)
		frare <- sapply(1:length(counts.sites),FUN=function(x) {round(table(unlist(strsplit(rep(names(counts.sites[[x]]),counts.sites[[x]]),split="/")))["1"]/sum(table(unlist(strsplit(rep(names(counts.sites[[x]]),counts.sites[[x]]),split="/")))),digits=4)})
		genind.obj <- suppressWarnings(vcfR::vcfR2genind(vcf.obj))
		if(is.na(popmap.path)){
			pops <- data.frame(V1=colnames(vcf.obj@gt[,-1]),V2=1)
		} else {
			pops <- read.table(popfile, header=F)
		}
		genind.obj@pop <- as.factor(pops$V2[match(colnames(vcf.obj@gt[,-1]),pops[,1])])
		# length(table(gt))
		# tab.df   <- as.data.frame(cbind(pop=genind.obj$pop,genind.obj@tab))
		# stats.wc  <- hierfstat::wc(tab.df)
		bs.nc   <- hierfstat::basic.stats(genind.obj)
		res.temp  <- data.frame(params=basename(mMNndirs[[i]]),as.data.frame(t(bs.nc$overall)))
		# res.temp$Fst.wc <- 
		res.temp$max.missingness.indv <- round(max(missingness.indv),digits=4)
		res.temp$min.missingness.indv <- round(min(missingness.indv),digits=4)
		res.temp$mean.missingness.indv <- round(mean(missingness.indv),digits=4)
		res.temp$max.missingness.site <- round(max(missingness.sites),digits=4)
		res.temp$min.missingness.site <- round(min(missingness.sites),digits=4)
		res.temp$mean.missingness.site <- round(mean(missingness.sites),digits=4)
		res.temp$frare.mean <- mean(frare)
		res.temp$nsite  <- nrow(gt)
		res.temp$nloci  <- length(unique(vcf.obj@fix[,"CHROM"]))
		res <- rbind(res,res.temp)
	}
	params.mat1  <- do.call(rbind,strsplit(res$params,"_"))
	params.mat2  <- gsub("^.","", params.mat1)
	mode(params.mat2) <- "integer"
	colnames(params.mat2) <- c("m","M","N","n")
	params.df <- as.data.frame(params.mat2)
	res2 <- data.frame(res,params.df)
	colnames(res2) <- c("params","Heterozygosity.observed","Heterozygosity.expected","Overall.GeneDiversity","Dst","Htp","Dstp","Fst","Fstp","Fis","Dest","max.missingness.indv","min.missingness.indv","mean.missingness.indv","max.missingness.site","min.missingness.site","mean.missingness.site","frare.mean","nsite","nloci","m","M","N","n")
	#rownames(res) <- basename(mMNndirs)
	#pshfiles <- list.files(dir,full.names=T,pattern="^populations[0-9]+.sh$")
	#file.remove(pshfiles)
	if(save){
		write.table(x=res2,file=file.path(save.in,"stacks_stats.txt"),quote=F,col.names=T,row.names=F,sep="\t")
	}
	res
}

#' @title summarize_VCFs
#' 
#' Same as 'summarize_stacks', except the full paths to input VCFs and an output directory are supplied as arguments.
#' 
#' @param VCF.paths Character string vector to one or more input VCF files.
#' @param save.as Path to location where output stats table should be saved. If NULL, the the table is not saved.
#' @param popmap.path Character string to population map. If NULL, all individuals in are assigned to the same 'population'.
#' @param include.out Character string vector indicating which output files to generate.
#' @return data frame(s) with basic stats
#' @export summarize_VCFs
summarize_VCFs <- function(VCF.paths,save.as,popmap.path=NULL,include.out=c("total","persite")){
	if(is.null(save.as)){
		save <- FALSE
	} else {
		save <- TRUE
	}
	if(is.null(popmap.path)){
		popmap.path <- NA
	}
	if("total" %in% include.out){
		res <- data.frame()
	} else {
		res <- NULL
	}
	if("persite" %in% include.out){
		res.sites <- data.frame()
	} else {
		res.sites <- NULL
	}
	for(i in 1:length(VCF.paths)){
		vcf <- VCF.paths[i]
		message(sprintf("processing '%s' [%s/%s]",basename(vcf),i,length(VCF.paths)))
		if(!file.exists(vcf)){
			stop(sprintf("could not find '%s'",vcf))
		}
		vcf.obj   <- vcfR::read.vcfR(vcf,verbose=F,checkFile=F)
		gt        <- gsub(":.+","",vcf.obj@gt[,-1])
		fx        <- vcf.obj@fix
		#gt.counts <- table(gt,useNA='always')
		#counts.sites <- sapply(X=1:nrow(gt), FUN=function(x) {A=table(gt[x,], useNA='always'); B=as.numeric(A[is.na(names(A))])/sum(A); B})
		counts.sites <- lapply(X=1:nrow(gt), FUN=function(x) {table(gt[x,], useNA='always')})
		#counts.indv <- sapply(X=1:ncol(gt), FUN=function(x) {A=table(gt[,x], useNA='always'); B=as.numeric(A[is.na(names(A))])/sum(A); B})
		counts.indv <- lapply(X=1:ncol(gt), FUN=function(x) {table(gt[,x], useNA='always')})
		# Fraction of individuals missing at each site
		missingness.sites <- round(sapply(X=1:length(counts.sites),FUN=function(x) {A=counts.sites[[x]]; as.numeric(sum(A[names(A) %in% c("./.",NA)]))/sum(A)}),digits=4)
		# Fraction of sites missing for each individual
		missingness.indv  <- round(sapply(X=1:length(counts.indv),FUN=function(x) {A=counts.indv[[x]]; as.numeric(sum(A[names(A) %in% c("./.",NA)]))/sum(A)}),digits=4)
		counts.sites.nonmissing <- lapply(X=1:length(counts.sites), FUN=function(x) {counts.sites[[x]][names(counts.sites[[x]]) %in% c("0/0","0/1","1/1")]})
		### From Nei 1987
		# Observed heterozygosity (= sHo = mHo in the basic.stats function)
		Ho.sites <- sapply(X=1:length(counts.sites),FUN=function(x){A=counts.sites[[x]]; B=sum(A[names(A) %in% c("0/0","1/1")])/sum(A[names(A) %in% c("0/0","0/1","1/1")]); 1-B})
		#Ho.total <- 1-mean(sapply(X=1:length(counts.sites),FUN=function(x){A=counts.sites[[x]]; sum(A[names(A) %in% c("0/0","1/1")])/sum(A[names(A) %in% c("0/0","0/1","1/1")])}))
		Ho.total <- mean(Ho.sites,na.rm=T)
		alleles.counts.sites   <- do.call(rbind,lapply(X=1:length(counts.sites.nonmissing),FUN=function(x){A <- counts.sites.nonmissing[[x]]; table(unlist(strsplit(rep(names(A),A),"/")))}))
		# x2i.bar  <- lapply(X=1:nrow(alleles.counts.sites),FUN=function(x){(alleles.counts.sites[x,]*alleles.counts.sites[x,])/sum(alleles.counts.sites[x,])})
		### xi = 'p' in the function basic.stats
		xi       <- do.call(rbind,lapply(X=1:nrow(alleles.counts.sites),FUN=function(x){(alleles.counts.sites[x,])/sum(alleles.counts.sites[x,])}))
		x2i.bar  <- xi*xi
		# x2.bar = 'sp2' = msp2 = mp2 in the function basic.stats
		x2.bar   <- apply(x2i.bar,1,sum)
		#x2i.bar  <- lapply(X=1:nrow(alleles.counts.sites),FUN=function(x){(alleles.counts.sites[x,])/sum(alleles.counts.sites[x,])*(alleles.counts.sites[x,])/sum(alleles.counts.sites[x,])})
		nk.sites <- sapply(1:length(Ho.sites),FUN=function(x){sum(alleles.counts.sites[x,])})
		# Same as 'mn' in basic.stats
		n <- nk.sites/2
		# Hs <- n/(n - 1) * (1 - sp2 - sHo/2/n)
		# Expected heterozygosity (Hs) at each site (= 'Hs' in basic.stats function)
		Hs.sites   <- n/(n-1)*(1-x2.bar-Ho.sites/2/n)
		# Expected heterozygosity (overall)
		Hs.total   <- mean(Hs.sites, na.rm=T)
		# Inbreeding coefficient at each site
		Fis.sites  <- 1-round((Ho.sites/Hs.sites),5)
		# Inbreeding coefficient (total)
		Fis.total  <- 1-round((Ho.total/Hs.total),5)
		mHs        <- n/(n-1)*(1-x2.bar-Ho.sites/2/n)
		Ht.sites   <- (1 - x2.bar + (mHs/n) - (Ho.sites/2/n))
		#Ht.total   <- (1 - x2.bar + (mHs/n) - (Ho.sites/2/n))
		Ht.total <- mean(Ht.sites,na.rm=T)
		#mFis <- 1 - Ho.sites/mHs
		#Dst  <- Ht - mHs
		#Dstp <- np/(np - 1) * Dst
		frare <- sapply(1:length(counts.sites),FUN=function(x) {round(table(unlist(strsplit(rep(names(counts.sites[[x]]),counts.sites[[x]]),split="/")))["1"]/sum(table(unlist(strsplit(rep(names(counts.sites[[x]]),counts.sites[[x]]),split="/")))),digits=4)})
		#if(FALSE){
		#	genind.obj <- suppressWarnings(vcfR::vcfR2genind(vcf.obj))
		#	if(is.na(popmap.path)){
		#		pops <- data.frame(V1=colnames(vcf.obj@gt[,-1]),V2=1)
		#	} else {
		#		pops <- read.table(popfile, header=F)
		#	}
		#	genind.obj@pop <- as.factor(pops$V2[match(colnames(vcf.obj@gt[,-1]),pops[,1])])
		#	# length(table(gt))
		#	# tab.df   <- as.data.frame(cbind(pop=genind.obj$pop,genind.obj@tab))
		#	# stats.wc  <- hierfstat::wc(tab.df)
		#	bs.nc   <- hierfstat::basic.stats(genind.obj)
		#	res.temp  <- data.frame(filename=basename(vcf),as.data.frame(t(bs.nc$overall)))
		#}
		#res.temp.sites     <- data.frame(filename=basename(vcf),locus=fx[,"CHROM"],position=fx[,"POS"],proportion.indv.missing=missingness.sites,indv.missing=counts.sites.nonmissing,Ho=Ho.sites,Hs=Hs.sites,Fis=Fis.sites,frequency.rare=frare)
		if("persite" %in% include.out){
			res.temp.sites.mat <- cbind(fx[,c("CHROM","POS")], missingness.sites, sapply(counts.sites.nonmissing,sum), round(Ho.sites,digits=5), round(Hs.sites,digits=5), round(Ht.sites,digits=5), round(Fis.sites, digits=5), round(frare,digits=5))
			colnames(res.temp.sites.mat) <- c("locus","position","proportion.indv.missing","indv.sampled","Ho","Hs","Ht","Fis","frequency.rare")
			res.temp.sites.df  <- data.frame(filename=basename(vcf),as.data.frame(res.temp.sites.mat))
			mode(res.temp.sites.df$proportion.indv.missing) <- mode(res.temp.sites.df$Ho) <- mode(res.temp.sites.df$Hs) <- mode(res.temp.sites.df$Ht) <- mode(res.temp.sites.df$Fis) <- mode(res.temp.sites.df$frequency.rare) <- "numeric"
			mode(res.temp.sites.df$indv.sampled) <- "integer"
			#min.ind.df.list <- lapply(X=sort(unique(res.temp.sites.df$indv.sampled)),FUN=function(x) {res.temp.sites.df[which(res.temp.sites.df$indv.sampled >= x),]})
			#Fis.minIndv     <- sapply(X=1:length(min.ind.df.list), FUN=function(x) {1-(mean(min.ind.df.list[[x]]$Ho,na.rm=T)/mean(min.ind.df.list[[x]]$Hs,na.rm=T))})
			#names(Fis.minIndv) <- as.character(sort(unique(res.temp.sites.df$indv.sampled)))
			res.sites <- rbind(res.sites,res.temp.sites.df)
		}
		if("total" %in% include.out){
			res.temp           <- data.frame(filename=basename(vcf))
			# res.temp$Fst.wc <- 
			res.temp$Heterozygosity.observed <- Ho.total
			res.temp$Heterozygosity.expected <- Hs.total
			res.temp$Overall.GeneDiversity   <- Ht.total
			res.temp$Fis                     <- Fis.total
			res.temp$max.missingness.indv    <- round(max(missingness.indv,na.rm=T),digits=4)
			res.temp$min.missingness.indv    <- round(min(missingness.indv,na.rm=T),digits=4)
			#res.temp$mean.missingness.indv <- round(mean(missingness.indv),digits=4)
			res.temp$max.missingness.site    <- round(max(missingness.sites,na.rm=T),digits=4)
			res.temp$min.missingness.site    <- round(min(missingness.sites,na.rm=T),digits=4)
			res.temp$missingness.total       <- round(mean(missingness.sites,na.rm=T),digits=4)
			res.temp$frare.mean <- mean(frare,na.rm=T)
			res.temp$nsite  <- nrow(gt)
			res.temp$nloci  <- length(unique(vcf.obj@fix[,"CHROM"]))
			res <- rbind(res,res.temp)
		}
	}
	persite.df <- function(df){
		df.list               <- lapply(X=1:length(unique(df$filename)),FUN=function(x){A= df[which(df$filename==unique(df$filename)[x]),]; rownames(A) <- 1:nrow(A); B=data.frame(A,row.id=as.numeric(rownames(A))); B})
		snps.first            <- lapply(1:length(df.list),function(x) 1:nrow(df.list[[x]]) %in% match(unique(df.list[[x]]$locus),df.list[[x]]$locus))
		df.list2              <- lapply(1:length(df.list),function(z) {df.list[[z]][order(df.list[[z]]$locus, df.list[[z]]$indv.sampled, decreasing=c(FALSE,TRUE)),]})
		snps.best             <- lapply(1:length(df.list2),function(x){1:nrow(df.list2[[x]]) %in% df.list2[[x]][match(unique(df.list2[[x]]$locus),df.list2[[x]]$locus),"row.id"]})
		df.list3              <- lapply(1:length(df.list),function(x){data.frame(df.list[[x]],is.firstSNP=snps.first[[x]],is.bestSNP=snps.best[[x]])})
		df.temp2              <- do.call(rbind,df.list3)
		df.temp2$filename     <- factor(df.temp2$filename)
		df.temp2$indv.sampled <- factor(df.temp2$indv.sampled)
		params.list    <- strsplit(levels(df.temp2$filename),"_")
		params.strings <- sapply(1:length(params.list),function(x) {paste(params.list[[x]][grep("^[m,M,N,n][0-9]+$",params.list[[x]])],collapse=" ")})
		each.x <- as.numeric(table(df.temp2$filename))
		stacks.params  <- unlist(sapply(1:length(params.strings),function(x){rep(params.strings[x],each.x[x])}))
		res.df <- data.frame(df.temp2,stacks.params=stacks.params)
		res.df$stacks.params <- factor(res.df$stacks.params)
		res.df$indv.sampled  <- factor(res.df$indv.sampled)
		res.df
	}
	if("persite" %in% include.out){
		res.sites2 <- persite.df(res.sites)
	}
	#colnames(res) <- c("filename","Heterozygosity.observed","Heterozygosity.expected","Overall.GeneDiversity","Dst","Htp","Dstp","Fst","Fstp","Fis","Dest","max.missingness.indv","min.missingness.indv","max.missingness.site","min.missingness.site","missingness.total","frare.mean","nsite","nloci")
	if(save){
		if("total" %in% include.out){
			write.table(x=res,file=save.as, quote=F,col.names=T,row.names=F,sep="\t")
		}
		if("persite" %in% include.out){
			save.as.persite <- paste0(tools::file_path_sans_ext(save.as),"_persite.txt")
			write.table(x=res.sites2,file=save.as.persite, quote=F,col.names=T,row.names=F,sep="\t")
		}
	}
	list(res,res.sites2)
}

### #' @title batch setup process_radtags
### #' 
### #' For each set pool of prepocessed samples this function creates a bash job file to run process_radtags.
### #' The output of this function is a character string vector with paths to bash job files, which can be piped to the function 'rsbatch' to submit the jobs.
### #' 
### #' @param save.in Character string with directory where output directories and files should be 
### #' @param barcodes.path Character string with path to file containing barcodes
### #' @param file Character string to input file, if processing single-end reads or stiched paired-end reads
### #' @param f1 Character string with path to first of the two input files, if processing paired-end reads. Default NULL.
### #' @param f2 Character string with path to second of the two input files, if processing paired-end reads. Default NULL.
### #' @param filetype Input file type
### #' @param clean Whether or not to clean the data, removing any read with an uncalled base. Default FALSE.
### #' @param quality Whether or not to discard reads with low quality scores. Default FALSE.
### #' @param rescue Whether or not to rescue barcodes and RAD-tags. Default FALSE.
### #' @param truncate Truncate final read length. Default NULL (no truncation).
### ##' @param t Number of threads to use for job parallization. Default 16.
### #' @param barcode.options Character string vector with one or more of the following: "inline_null","index_null","inline_inline","index_index","inline_index","index_inline".
### #' @param process_radtags.path Character string with path to the process_radtags executable. Default NULL, in which case the path is derived from the default setting, which can be set with config_miscwrappers('<path/to/process_radtags>').
### #' @return Character string vector with the paths to bash job files.
### #' @export process_radtags_setup
### process_radtags_setup <- function(save.in,barcodes.path,file=NULL,f1=NULL,f2=NULL,filetype="gzfastq", paired=F,interleaved=F,barcode.options=c("inline_null"),clean=F,quality=F,rescue=F,truncate=NULL,D=F,E=NULL,w=0.15,s=10,y=NULL,process_radtags.path=NULL){
### 	if(is.null(process_radtags.path)){
### 		if(any(config_miscwrappers()[,"program"]=="process_radtags")){
### 			process_radtags.path <- config_miscwrappers()[config_miscwrappers()[,"program"]=="process_radtags","exe_path"]
### 		} else {
### 			stop("Set 'process_radtags.path' argument, or use config_miscwrappers(exe.path=<path/to/process_radtags>) to set the package default path for process_radtags")
### 		}
### 	}
### 	mMNndirs <- list.dirs(save.in,recursive=F)
### 	outdirs <- mMNndirs
### 	if(is.null(popmap)){
### 		if("popmap.txt" %in% list.files(save.in)){
### 			popmap.path <- file.path(save.in,"popmap.txt")
### 		} else {
### 			popmap.df <- data.frame(indv=basename(unique(gsub(".snps.tsv.gz","",list.files(mMNndirs[1],full.names=T,pattern=".snps.tsv.gz")))),pop=1)
### 			popmap.df2 <- popmap.df[popmap.df[,"indv"]!="catalog",]
### 			popmap.path <- file.path(save.in,"popmap.txt")
### 			write.table(popmap.df2,file=popmap.path,col.names=F,row.names=F,sep="\t",quote=F)
### 		}
### 	} else {
### 		popmap.path <- popmap
### 	}
### 	# process_radtags -P ./stacks/ -M ./popmap -R ./samplesdir -s sample1.fq.gz -s sample2.fq.gz-t 8
### 	if(!is.null){
### 		sample.paths <- paste("-s",samples)
### 	} else {
### 		sample.paths <- NA
### 	}
### 	if(!is.null){
### 		samples.dir <- paste("-R",samples.dir)
### 	} else {
### 		samples.dir <- NA
### 	}
### 	all.df   <- data.frame(t2b=sq(process_radtags.path),P=paste("-P",sq(outdirs)),M=paste("-M",popmap.path),R=samples.dir,s=sample.paths,p=paste("-t",t))
### 	stacks.strings <- sapply(1:nrow(all.df),FUN=function(x) {paste(all.df[x,],collapse=" ")})
### 	sh.lines  <- lapply(1:length(stacks.strings),FUN=function(x){c("#!/bin/bash",stacks.strings[x])})
### 	sh.paths  <- file.path(save.in,paste0("process_radtags",1:length(stacks.strings),".sh"))
### 	sh.write  <- lapply(1:length(sh.lines),FUN=function(x){writeLines(text=sh.lines[[x]],con=sh.paths[x])})
### 	sh.paths
### }
### 
### /panfs/pfs.local/work/bi/bin/stacks-2.5/process_radtags -i gzfastq -f /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R1_R2.stitched_v3.fq.gz --inline_null -b /panfs/pfs.local/home/j926w878/work/ddRAD/files-to-run_Stacks_pipeline/KTJW_batch4.barcodes.txt -o /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/process_radtags_output/ --renz_1 sbfI --renz_2 mspI -c -q -r -D -w 0.15 -s 20 --barcode_dist_1 2 2>/panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/process_radtags_output/KTJW.demux.out
####### Preprocess sequences
### Copies each forward (read1) and reverse (read2) files to files named "R1.fq.gz" and "R2.fq.gz" 
### cp /panfs/pfs.local/home/j926w878/work/ddRAD/raw_data/KTJW1/KT_JW_P2_4_CKDL200146889-1a-8_H7GWYBBXX_L2_1.fq.gz /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R1.fq.gz &
### cp /panfs/pfs.local/home/j926w878/work/ddRAD/raw_data/KTJW1/KT_JW_P2_4_CKDL200146889-1a-8_H7GWYBBXX_L2_2.fq.gz /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R2.fq.gz
### wait
### # unzips the copied files
### cd /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/; gunzip /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R1.fq.gz &
### cd /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/; gunzip /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R2.fq.gz
### wait
### ### for each R1 file, cut the last 2 bp (probably shitty quality)
### ### this assumes that each read is 150bp (usually true) so need to check that first
### awk 'NR%4==2' /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R1.fq | awk '{ print "\n\n\n"$1;}' | sed -e '1,2d' | cut -c 1-148 > /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R1.seq
### wait
### ### For each R1 file, cut the last 2 characters from the quality line (aka keep first 148 chars)
### awk 'NR%4==0' /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R1.fq| awk '{ print "\n\n\n"$1;}' | cut -c 1-148 | rev > /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R1.qual
### wait
### ### for each R2 file, cut the last 2 bp (probably shitty quality), and reverse complement it
### ### this assumes that each read is 150bp, so need to check that first
### awk 'NR%4==2' /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R2.fq | awk '{ print "\n\n\n"$1;}' | sed -e '1,2d' | cut -c 1-148 | rev | tr ACGT TGCA > /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R2.seq.rc
### wait
### ### for each R2 file, cut the last 2 characters from the quality line, and then reverse the order of the remaining qual scores to match the order of the reverse complemented DNA sequences
### awk 'NR%4==0' /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R2.fq| awk '{ print "\n\n\n"$1;}' | cut -c 1-148 | rev > /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R2.qual
### wait
### ## gets every-other line of R1, starting with first line, and pastes it onto the same line of a new file
### ## no need to run this for R2, because the description lines are the same in R1 and R2 files (and will be the same in the stitched reads file)
### awk 'NR%2==1' /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R1.fq | awk '{ print $0"\n";}' > /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R1_R2.descriptions
### ### stitch trimmed R1 to the trimmed rev. comp. of R2 to form a new "stitched" read
### ### includes the trimmed version of R1 (i.e., with last two base removed)
### ### does NOT include the 450nt insert region
### ######
### paste /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R1_R2.descriptions /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R1.seq /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R1.qual /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R2.seq.rc /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R2.qual | awk '{gsub("\t","",$0); print $0;}' | awk '{gsub(" 1:N:0:"," stitched:N:0:",$0); print $0;}' > /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R1_R2.stitched_v3.fq
### wait
### ### to save disk space, gzip the resulting stitched reads (STACKS can read gzipped files)
### ### does the same thing as gzip.sh, but for the files R1_R2.stitched_v3.fq
### cd /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/; gzip -9 /panfs/pfs.local/home/j926w878/scratch/ddRAD/raw_data/KTJW1/R1_R2.stitched_v3.fq
### wait


