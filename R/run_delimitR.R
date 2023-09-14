#' @title Run delimitR
#' 
#' Function to use delimitR to setup standard demographic models and run fastsimcoal2. Custom models other than the default ones constructed by delimitR can also be included.
#' 
#' @param fsc.path Path to fastsimcoal executable, or the name of the executable if it is on your path or in your current directory. Default is "fsc26", but this will be different if you have a different version of fastsimcoal.
#' @param working.dir Path to working directory where observed SFS (*.obs) and traits files (.
#' @param observedSFS Name of the observed SFS file (with or without extension). Default NULL.
#' @param traitsfile Name of the traits file (including extension). If 'observedSFS' is NULL and 'vcf' is non-NULL, the traits file must be formatted as required by easySFS; otherwise, the traits file must be formatted like as required by fastsimcoal2.
#' @param easySFS.path Optional character string with path to easySFS executable. Only required when 'observedSFS' is NULL and 'vcf' is non-NULL.
#' @param vcf Optional character string with path to VCF file. Default is NULL. If this is set, the value of observedSFS will be ignored and the an SFS will be created by calling easySFS. The format of the traits file must be as required by easySFS; the fsc2-format traits file will automatically be generated to pass to delimitR and fsc2.
#' @param custom_models Character string with path or name (if in current directory) to a directory containing custom models to use in addition to the standard models constructed by delimitR. Default is to not include any custom models.
#' @param model.subset NULL (default) or a numerical vector indicating which models to include. If NULL, all models generated are used. If a numerical vector, each number corresponds to the number after the prefix and before the extension of each .tpl and .est file.
#' @param obsprefix Character string or an integer in c(0, 1). Character string specifies the prefix to use for .est and .tpl files defining demographic models. These files must have the form "obsprefix_n.tpl" or "obsprefix_n.est", where n is an integer in 1:N models.
#' If 0 (the default), "model" will be used as the prefix.
#' If 1, the function will attempt to infer the prefix from filenames of custom models supplied. Do not set to 1 if you are not supplying custom models. This will also fail if multiple prefixes are used for custom models, or if .est and .tpl files do not have the form 'obsprefix_n.est' and 'obsprefix_n.tpl'.
#' @param observedtree.manyspecies Character string describing the observed species tree in Newick format. Use NULL (the default) unless > 3 populations (species) are observed.
#' @param divwgeneflow If models should be created that include divergence with gene flow. Default is TRUE.
#' @param seccontact If models should be created that include secondary contact. Default is TRUE.
#' @param maxedges Number of migration events to include in standard delimitR models. At present must equal 1. You can create custom models if you want to include multiple edges (migration events).
#' @param obspopsizeprior List of numeric vectors indicating the min and max prior values for each population. Default is "default", which uses 10000 and 100000 for the minimum and maximum, respectively, for all populations. You probably shouldnt use "default".
#' The input has the form: list(c(min_Pop0,max_Pop0),c(min_Pop1,max_Pop1), ...)
#' Note that these are in terms of the number of haploid individuals, so use n*2 for diploid species.
#' @param obsdivtimeprior List of numeric vectors indicating min and max values for divergence times (units = number of generations). Default is for a two species case (single divergence), and the values should be changed given your situation.
#' Divergence time priors should be provided in order of coalescent interval.
#' @param obsmigrateprior List containing one numeric vector with min and max values of migration rates in the migration matrix. A single non-zero migration matrix is allowed in the standard delimitR models. This prior does not affect any custom models in the custom_models directory.
#' @param nreps Number of fastsimcoal replicates. Default is 10000.
#' @param nRFtrees Number of RF trees to use. Default is 500.
#' @param nclasses.bin Number of bins for turning mSFS into binned bSFS. This number should not be greater than the sample size of the population with the fewest samples, as this results in sparse sampling of the SFS. Default is 5.
#' Large values lead to a more complete summary of the data, but also lead to a more sparsely sampled SFS and increased computation times.
#' @param ncores Number of cores to use. Default is 2.
#' @param run.parallel Whether or not simulations of different models should be run in parallel. Default is TRUE, but in some cases weird issues arrise and in it can be useful to change to FALSE.
#' @param threshold Number indicating threshold used for downsampling during SFS construction if this was done. Default is 100 (no downsampling).
#' @param clean.wd Whether or not to delete extra/temporary files after running. Default is TRUE. Changing to FALSE can help with debugging but you will regret doing setting to FALSE for long-term because many files are generated.
#' @param shuffle.traits Whether or not to randomly shuffle the assignment of traits to individuals. Default is FALSE. This option is not yet implemented but will probably require calling easySFS, which will require data input as vcf and some way to automate SFS options.
#' @return NULL. Output files are written working.dir
#' @export run_delimitR
run_delimitR <- function(fsc.path="fsc26", working.dir=getwd(), observedSFS=NULL,traitsfile, easySFS.path=NULL, vcf=NULL, custom_models=NULL, model.subset=NULL, obsprefix=0, observedtree.manyspecies=NULL, divwgeneflow = TRUE, seccontact = TRUE, obsdivtimeprior = list(c(50000,100000)), obsmigrateprior = list(c(0.000005,0.00005)), obspopsizeprior = "default", maxedges = 1, nreps=10000,nRFtrees=500,nclasses.bin=5,ncores=2,run.parallel=T,threshold=100,clean.wd=T,shuffle.traits=F){
	# Create MSFS and fastsimcoal2-style traits file if 'vcf' and 'easySFS.path' are non-NULL & 'traitsfile' holds path to easySFS-style traits file & observedSFS is NULL
	if(is.null(observedSFS) & !is.null(vcf) & !is.null(easySFS.path)){
		data.paths   <- easySFS_fsc2(easySFS.path=easySFS.path,vcf=vcf,traits=traitsfile,output.dir=working.dir)
		observedSFS  <- basename(data.paths["MSFS"])
		traitsfile   <- basename(data.paths["traits"])
	}

	# set working directory. Not sure if this is good to do within a function.
	# setwd(working.dir)
	# use basename and remove extension
	observedSFS <- tools::file_path_sans_ext(basename(observedSFS))
	# sample sizes inferred from line 2 of the observedSFS
	sampleSizes   <- data.table::fread(paste0(observedSFS,".obs"),skip=1,nrows=1)
	# snps spectrum from line 3 of the observedSFS
	snps.spectrum <- data.table::fread(paste0(observedSFS,".obs"),skip=2)
	# set the prefix used for models
	if(obsprefix==0){
		obsprefix <- "model"
	} else {
		if(obsprefix==1){
			# There are probably other ways that this could go wrong.
			if(is.null(custom_models)){
				stop("custom_models must be non-NULL if obsprefix is 1")
			} else {
				tpl.custom <- list.files(custom_models,pattern = ".*\\.tpl")
				est.custom <- list.files(custom_models,pattern = ".*\\.est")
				tplest.custom <- c(tpl.custom,est.custom)
				prefix.custom <- gsub("_.+","",tplest.custom)
				if(length(unique(prefix.custom))>1){
					stop("models in custom_models directory must have the same prefix")
				} else {
					obsprefix <- prefix.custom
				}
			}
		}
	}
	# Number of populations (or species).
	obsspecies  <- as.numeric(sampleSizes[[1]])
	
	if(obsspecies %in% 1:3){
		observedtree <- c('(0);','(0,1);','((0,1),2);')[obsspecies]
	} else {
		if(is.null(observedtree.manyspecies)){
			stop("observedtree.manyspecies must be non-NULL when >3 populations (species) are defined in traits file")
		} else {
			observedtree <- observedtree.manyspecies
		}
	}
	
	### Creating a symmetric migration matrix used in standard delimitR models
	migmatrix       <- matrix(data=TRUE,nrow=obsspecies,ncol=obsspecies)
	diag(migmatrix) <- FALSE
	
	## The number of "alleles" retained after downsampling mSFS
	# These are the values on line 2 of mSFS file except for the first value
	obssamplesize  <- as.numeric(unlist(strsplit(sampleSizes[[2]],split=" ")))
	
	## Number of linkage blocks to simulate.
	# For unlinked SNPs, this is equal to the number of SNPs used to build your observed SFS, which is the sum of the values on line 3 of the mSFS file!
	obssnps <- sum(snps.spectrum)
	# priors on population sizes.
	if(!is(obspopsizeprior,"list")){
		if(obspopsizeprior=="default"){
			obspopsizeprior <- rep(list(c(10000,100000)), obsspecies)
		} else {
			stop("obspopsizeprior must be a list of numerical vectors or the character string 'default'")
		}
	}
	## Set up your prior models for fastsimcoal2
	## Produces standard set of delimitR models. These are held in one .tpl and one .est file for each model. Usually four models if two populations assigned in traits file.
	delimitR::setup_fsc2(tree=observedtree,nspec=obsspecies,samplesizes=obssamplesize,nsnps=obssnps,prefix=obsprefix,migmatrix=migmatrix,popsizeprior=obspopsizeprior,divtimeprior=obsdivtimeprior,migrateprior=obsmigrateprior,secondarycontact= seccontact,divwgeneflow= divwgeneflow,maxmigrations = maxedges)
	## Copies custom models (if supplied) into working directory and checks/edits 'sample sizes' and 'number of independent loci' to match obs file.
	if(!is.null(custom_models)){
		# paths to custom models
		old.paths  <- list.files(custom_models,full.names=T)
		# New names to use for custom_models to be consistent with the prefix used for the standard set of models
		new.names  <- gsub("^.+_",paste0(obsprefix,"_"),basename(old.paths))
		# Full path where custom models will be copied to
		new.paths  <- file.path(working.dir,new.names)
		# Copy the custom models to the working directory and use the new names.
		file.copy(from=old.paths,to=new.paths)
		#tpl.files <- list.files(pattern=".tpl$")
		tpl.files <- new.paths[grep(".tpl$",new.paths)]
		for(i in 1:length(tpl.files)){
			tpl.lines.i <- suppressWarnings(readLines(tpl.files[i]))
			tpl.lines.i[(match("//Sample sizes",tpl.lines.i)+1):(match("//Sample sizes",tpl.lines.i)+length(obssamplesize))] <- obssamplesize
			tpl.lines.i[match("//Number of independent loci [chromosome]",tpl.lines.i)+1] <- obssnps
			writeLines(text=tpl.lines.i,con=tpl.files[i])
		}
	}
	### Only consider models in model.subset
	if(!is.null(model.subset)){
		# filenames of all models in working directory
		tpl.all <- list.files(working.dir,pattern = ".*\\.tpl$")
		est.all <- list.files(working.dir,pattern = ".*\\.est$")
		# filenames of all models to consider
		tpl.use <- paste0(obsprefix,"_",model.subset,".tpl")
		est.use <- paste0(obsprefix,"_",model.subset,".est")
		# filepaths of models in current directory that should not be considered
		files.remove  <- file.path(working.dir,setdiff(c(tpl.all,est.all),c(tpl.use,est.use)))
		# removing model files for the models to not use
		system(paste("rm",files.remove,collapse=" & "))
		# Rename models such that they are numbered sequentially from 1, and make a table that shows oldname vs. newname
		new.tpl.extensions <- paste0("_",1:length(tpl.use),".tpl")
		new.est.extensions <- paste0("_",1:length(tpl.use),".est")
		tpl.newnames       <- sapply(X=1:length(tpl.use),FUN=function(x){gsub("_.*\\.tpl$",new.tpl.extensions[x],tpl.use[x])})
		est.newnames       <- sapply(X=1:length(est.use),FUN=function(x){gsub("_.*\\.est$",new.est.extensions[x],est.use[x])})
		system(paste0(paste("mv",file.path(working.dir,tpl.use),file.path(working.dir,tpl.newnames)),collapse=" & "))
		system(paste0(paste("mv",file.path(working.dir,est.use),file.path(working.dir,est.newnames)),collapse=" & "))
		### make and then save a table showing oldnames and newnames of models
		modelnames_table  <- cbind(oldnames=gsub(".tpl","",tpl.use),newnames=gsub(".tpl","",tpl.newnames))
		write.table(modelnames_table,file="modelnames_table.txt",quote=F,sep="\t",row.names=F)
		# Info from the modelnames table will be used later to revert back to oldnames once the analysis is complete
	}
	# Run fastsimcoalsims in parallel
	fastsimcoalsims_Par(prefix=obsprefix,pathtofsc=fsc.path,nreps=nreps,ncores=ncores,run.parallel=run.parallel)
	# Specify bin number for turning mSFS into binned bSFS
	nclasses <- nclasses.bin
	# Make the full prior
	FullPrior <- delimitR::makeprior(prefix=obsprefix,nspec=obsspecies,nclasses=nclasses,mydir=getwd(),traitsfile = traitsfile,threshold=threshold, thefolder = "Prior",ncores = ncores)
	
	## Created the reduced prior by removing invariant data.
	ReducedPrior <- delimitR::Prior_reduced(FullPrior)
	# Switching back to the original model names if a subset of models were used
	if(!is.null(model.subset)){
		model.indices <- lapply(X=1:nrow(modelnames_table),FUN=function(x){which(FullPrior[,"Model"] == paste0(modelnames_table[x,"newnames"],"_MSFS"))})
		model.names   <- unlist(lapply(X=1:length(model.indices),FUN=function(x){rep(paste0(modelnames_table[x,"oldnames"],"_MSFS"),length(model.indices[[x]]))}))
		model.names   <- model.names[order(unlist(model.indices))]
		FullPrior$Model <- as.factor(model.names)
		### Also change filenames back to original names, and use these names for the binned files in the priors folder
		tpllist    <- system(paste0("ls ", obsprefix, "*.tpl"),intern = T)
		estlist    <- system(paste0("ls ", obsprefix, "*.est"),intern = T)
		different.names.table <- modelnames_table[modelnames_table[,"newnames"] != modelnames_table[,"oldnames"],,drop=F]
		system(paste(paste("mv",file.path(working.dir,paste0(different.names.table[,"newnames"],".tpl")), file.path(working.dir,paste0(different.names.table[,"oldnames"],".tpl"))),collapse=" & "))
		system(paste(paste("mv",file.path(working.dir,paste0(different.names.table[,"newnames"],".est")), file.path(working.dir,paste0(different.names.table[,"oldnames"],".est"))),collapse=" & "))
		system(paste(paste("mv",file.path(paste0(working.dir,"/Prior"),paste0("Binned_Processed_",different.names.table[,"newnames"],"_MSFS.obs")), file.path(paste0(working.dir,"/Prior"),paste0("Binned_Processed_",different.names.table[,"oldnames"],"_MSFS.obs"))),collapse=" & "))
	}

	## Build random forest generated from the full and reduced priors
	myRF <- delimitR::RF_build_abcrf(ReducedPrior=ReducedPrior,FullPrior=FullPrior,ntrees=nRFtrees)
	classification.error <- myRF[["model.rf"]]$confusion.matrix
	write.csv(classification.error, file = "classification.error.csv")
	pdf("myRF.pdf")
	plot(myRF, training = ReducedPrior)
	;
	dev.off()
	# Prep observed data
	myobserved <- delimitR::prepobserved(observedSFS,FullPrior,ReducedPrior,nclasses,obsspecies,traitsfile=traitsfile,threshold = threshold)
	# Apply Random Forest classifier to observed data
	prediction <- delimitR::RF_predict_abcrf(RF=myRF, Observed_Reduced=myobserved, ReducedPrior=ReducedPrior, FullPrior=FullPrior, ntrees=nRFtrees)
	#### Need to update names of headers in  the prediction matrix if you used a subset of models!
	## Also need to verify that the model names are correctly reassigned because sorting names is weird such that "model10" will be sorted earlier than "model2"
	prediction2 <- as.matrix(prediction)
	### This part doesn't quite work.
	if(!is.null(model.subset)){
		### may be as simple as this:
		colnames(prediction2) <- c("selected model", paste(sort(modelnames_table[,"oldnames"]),"votes"),"post.prob")
		# Then need to sort all but first and last columns such that models are sorted by integer value rather than by each character in name.
#		# selected model votes model1 votes model2 votes model3 votes model4 votes model5 votes model6 votes model7 post.prob
#		prediction2 <- prediction2[c(1,(order(as.numeric(gsub(paste0(obsprefix,"_"),"",sort(modelnames_table[,"oldnames"]))))+1),length(modelnames_table))]
	}
	# Save prediction. This table shows which model was best supported using the RF
	write.csv(prediction2, file="prediction.csv")
	
	# Read prediction table, and then rewrite 
	#
	### Save a copy of the model specifications
	# tpl and est filenames
	tpllist    <- system(paste0("ls ", obsprefix, "*.tpl"),intern = T)
	# model info stored in a character vector. Will fail if 'incomplete final lines'
	model.info <- unlist(lapply(X=1:length(tpllist),FUN=function(x){c(paste("model",x),readLines(tpllist[x]),"")}))
	## Write out model info to file
	writeLines(model.info,con="models.txt")
	# Remove extra/temporary files created during analysis.
	if(clean.wd){
		delimitR::clean_working(prefix=obsprefix)
	}
}
#' @examples
#' 
#' run_delimitR(vcf="PATH/TO/VCF",traitsfile="PATH/TO/TRAITSFILE")


#' @title fastsimcoalsims_Par
#' 
#' Function to run fastsimcoalsims in parallel. This is modified from Paul Hime's function of the same name. Implemented in the function rundelimitR.
#' Simulates the site frequency spectrum for a set of input models each defined by their .est and .tpl files.
#' 
#' @param prefix Prefix used for .tpl and .est model files
#' @param pathtofsc Path to the fastsimcoal2 executable
#' @param nreps Number of replicate runs to do for each model. Default is 10000.
#' @param wd Working directory. Default is to use the current working directory.
#' @param ncores Number of cores to use
#' @param run.parallel If simulations on different models should also be run in parallel. Default is TRUE. In some cases running in parallel causes weird issues and in such cases try changing this to FALSE.
#' @return NULL. Results are written to working directory.
#' @export fastsimcoalsims_Par
fastsimcoalsims_Par <- function(prefix, pathtofsc, nreps=10000,wd=getwd(),ncores=2,run.parallel=T){
	#tpllist     <- system(paste0("ls ", prefix, "*.tpl", sep = ""),intern = T)
	#estlist     <- system(paste0("ls ", prefix, "*.est", sep = ""),intern = T)
	tpllist      <- list.files(pattern = ".*\\.tpl$")
	estlist      <- list.files(pattern = ".*\\.est$")
	commands     <- unlist(lapply(X=1:length(tpllist),FUN=function(j){paste0(pathtofsc, " -t ", file.path(wd,paste0(prefix, "_", j, ".tpl"))," -e ", file.path(wd,paste0(prefix, "_", j, ".est")), " -n 1 --msfs --cores ",ncores," -q --multiSFS -x -E",nreps)}))
	commands.par <- paste(commands,collapse=" & ")
	if(run.parallel){
		system(commands.par, ignore.stdout = TRUE,wait=TRUE)
	} else {
		for(i in 1:length(commands)){
			system(commands[i], ignore.stdout = TRUE,wait=TRUE)
		}
	}
}


#' Run easySFS even easier to generate input MSFS and population files for fastsimcoal2
#'
#' Given VCF and traits files usually needed by easySFS, generate the fsc2 MSFS and the pops file needed for fastsimcoal2
#' Chooses the projection that yields the greatest number of segregating sites, without trying to balance for the number of individuals.
#' 
#' IMPORTANT: You must include the next two lines in the sh file for any job that uses this function. For now this function will only work for me.
#' module load anaconda/4.7
#' conda activate easySFS
#' 
#' requires gtools package
#' 
#' @param easySFS.path Character string with path to the easySFS executable
#' @param vcf Character string with path to the VCF file
#' @param traits Character string with path to the traits file
#' @param output.dir Character string with path to directory where the two output files should be saved.
#' @export easySFS_fsc2
easySFS_fsc2 <- function(easySFS.path,vcf,traits,output.dir){
	prev.tempfile  <- tempfile()
	prev.command   <- paste(easySFS.path,"-i",vcf,"-p",traits,"--preview >",prev.tempfile)
	run.prev       <- system(prev.command)
	prev.lines     <- readLines(prev.tempfile)
	traits.df      <- read.table(traits,header=F)
	popsIDs        <- sort(unique(traits.df[,2]))
	popheaderlines <- match(popsIDs,prev.lines)
	popprojlines   <- popheaderlines + 1
	proj.strings   <- prev.lines[popprojlines]
	num.indv       <- NULL
	for(i in 1:length(proj.strings)){
		split1.i <- unlist(strsplit(proj.strings[i],split="\t"))
		split1.i <- gsub("(","",split1.i,fixed=T)
		split1.i <- gsub(")","",split1.i,fixed=T)
		split1.i <- gsub(" ","",split1.i,fixed=T)
		split2.i <- lapply(X=strsplit(split1.i,split=",",fixed=T),FUN=as.numeric)
		mat.i    <- do.call(rbind,split2.i)
		ordered.mat.i <- mat.i[order(mat.i[,2],decreasing=T),]
		num.indv <- c(num.indv,ordered.mat.i[1,1])
	}
	names(num.indv) <- popsIDs
	num.indv.string <- paste(num.indv,collapse=",")
	num.indv.string2 <- paste(num.indv,collapse="_")
	#temp.outdir     <- file.path(tempdir(),paste(sample(letters,10),collapse=""))
	temp.outdir     <- file.path(output.dir,paste(sample(letters,10),collapse=""))
	easySFS.command <- paste(easySFS.path,"-i",vcf,"-p",traits,"--proj",num.indv.string,"-o",temp.outdir)
	system(easySFS.command)
	### temporary filepath where the MSFS file exists
	output1.temp.path <- list.files(file.path(temp.outdir,"fastsimcoal2"),pattern="_MSFS.obs$",full.names=TRUE)
	### Now time to make the  fastsimcoal2 style traits file
#	if(diploid){
		#traits.col1 <- unlist(lapply(X=1:length(num.indv),FUN=function(x){gtools::mixedsort(c(paste0("pop",LETTERS[x],1:num.indv[x],"_a"),paste0("pop",LETTERS[x],1:num.indv[x],"_b")))}))
#		traits.col1 <- c(unlist(lapply(X=1:length(num.indv),FUN=function(x){gtools::mixedsort(c(paste0("pop",popsIDs[x],"ind",1:num.indv[x],"_a"),paste0("pop",popsIDs[x],"ind",1:num.indv[x],"_b")))})))
#		traits.col2 <- c(unlist(sapply(X=1:length(num.indv),FUN=function(x){rep(popsIDs[x],(num.indv[x])*2)})))
#	} else {
	traits.col1 <- c(unlist(lapply(X=1:length(num.indv),FUN=function(x){gtools::mixedsort(paste0("pop",popsIDs[x],"ind",1:num.indv[x]))})))
	traits.col2 <- c(unlist(sapply(X=1:length(num.indv),FUN=function(x){rep(popsIDs[x],num.indv[x])})))
#	}
	traits.df <- data.frame(traits=traits.col1,species=traits.col2)
	### Names for output files
	projprev.filename <- "easySFS_projection_preview.txt"
	msfs.filename     <- paste0(tools::file_path_sans_ext(basename(vcf)),"_K",length(popsIDs),"_",num.indv.string2,"_MSFS.obs")
	traits.filename   <- paste0(tools::file_path_sans_ext(basename(vcf)),"_K",length(popsIDs),"_",num.indv.string2,"_TRAITS.txt")
	### Move easySFS projection preview and MSFS files, and write traits file to output.dir.
	mv.cmd1 <- paste("mv",output1.temp.path,file.path(output.dir,msfs.filename))
	mv.cmd2 <- paste("mv",prev.tempfile,file.path(output.dir,projprev.filename))
	system(mv.cmd1)
	system(mv.cmd2)
	write.table(traits.df,file.path(output.dir,traits.filename),col.names=T,row.names=F,sep="\t",quote=F)
	result.paths <- c(MSFS=file.path(output.dir,msfs.filename),traits=file.path(output.dir,traits.filename))
	### Remove temporary files
	if(dir.exists(temp.outdir)){
		system(paste("rm -R",temp.outdir))
	}
	### Output is a character vector holding the paths to location where MSFS and traits files were written
	result.paths
}

