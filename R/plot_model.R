#' @title Plot Demographic Model
#' 
#' This function creates a graphical representation of a demographic model defined by a template (.tpl) and estimation (.est) file.
#' Concern: Setting initial growthrate values to a variable name will probably lead to failure.
#' fastsimcoal features not implemented: 'nomig' historical event directive.
#' 
#' @param tpl.path Path to tpl file of model
#' @param est.path Path to est file of model
#' @param model.title Character string with the title text to use above the plot, e.g., the model name/ID. Default is "".
#' @param popnames Optional character vector with label names to use for populations. Default is NULL (don't print names on plot); supplying a single integer 1 will use names of the form "pop0", "pop1",... If a character vector is supplied, the number of names supplied should equal the number of populations indicated on line two of template file. Ordering of names correspond to population 0, population 1, ..., as used in the template and estimation files.
#' @param show.priors Whether or not to also plot the prior probability distributions of parameters. Default is FALSE. This is still in progress.
#' @param pops.spacer Number from 1 to 2 indicating amount of white space to allow between populations. Default is 1.05.
#' @param plotmargin Numerical vector indicating how much extra plotting area to add below, left, above, and rigt of the graph.
#' @param label.cex Size to use for text labels (not including axis labels)
#' @param warn.about.comments Whether or not to warn about in-line comments found in unusual places of '.tpl' file. Default is FALSE.
#' @param priors.panel.width Fraction of the plot to use for plotting the priors panel. Default is 0.3. Ignored if show.priors is FALSE. Setting this argument to zero is the same as setting show.priors FALSE.
#' @param use.mean Whether or not to use the mean or a random sample from each paramater's prior distribution when drawing graphs. Default is TRUE.
#' @return An object of class 'recordedplot' that stores a copy of the plot. See examples.
#' @export plot_model
plot_model <- function(tpl.path, est.path, model.title="", popnames=NULL, show.priors=TRUE,pops.spacer=1.1, plotmargin=c(1.1,1.01,1.05,1), label.cex=0.75, warn.about.comments=FALSE,priors.panel.width=0.3,use.mean =TRUE){
	##### Defaults for debugging
	if(!exists("model.title")){
		model.title <- ""
	}
	if(!exists("popnames")){
		popnames <- NULL
	}
	if(!exists("show.priors")){
		show.priors <- TRUE
	}
	if(!exists("pops.spacer")){
		pops.spacer <- 1.1
	}
	if(!exists("plotmargin")){
		plotmargin=c(1.1,1.01,1.05,1)
	}
	if(!exists("label.cex")){
		label.cex=0.75
	}
	if(!exists("warn.about.comments")){
		warn.about.comments=FALSE
	}
	if(!exists("priors.panel.width")){
		priors.panel.width=0.3
	}
	if(!exists("use.mean")){
		use.mean=TRUE
	}

	########
	tpl.lines       <- readLines(tpl.path, warn=F)
	est.lines.temp  <- readLines(est.path, warn=F)
	## drops comment lines and empty lines
	est.lines         <- est.lines.temp[-c(grep("^//",est.lines.temp),which(est.lines.temp==""))]
	rules.header.line <- grep("^.RULES",est.lines)
	### Lines in est file that start with an upper or lowercase letter define rules
	est.rules.lines     <- as.list(est.lines[grep("^[letters,LETTERS]",est.lines)])
	### Lines that start with 0 or 1 and occur earlier than the RULES header line define simple parameters
	est.parameters.lines <- as.list(est.lines[grep("^[0-1]",est.lines[1:rules.header.line])])
	### Lines that start with 0 or 1 and occur later than the RULES header line define complex parameters
	est.complex.lines   <- as.list(est.lines[grep("^[0-1]",est.lines[rules.header.line:length(est.lines)]) + (rules.header.line-1)])
	if(length(est.parameters.lines)>0){
		variable.names <- unlist(lapply(X=est.parameters.lines,FUN=function(x){unlist(strsplit(x,split=" +"))[2]}))
		if(length(est.complex.lines)>0){
			variable.names <- c(variable.names,unlist(lapply(X=est.complex.lines,FUN=function(x){unlist(strsplit(x,split=" +"))[2]})))
		}
	}
	if(length(est.rules.lines)>0){
		rules.exist <- TRUE
		est.rules.lines <- gsub(" ","",est.rules.lines)
		for(i in 1:length(variable.names)){
			est.rules.lines <- gsub(variable.names[i],paste0(" ",variable.names[i]," "),est.rules.lines,fixed=T)
		}
		est.rules.lines <- gsub("^ +","",est.rules.lines)
		est.rules.lines <- as.list(gsub(" +$","",est.rules.lines))
		rules.mat       <- do.call(rbind,lapply(X=est.rules.lines,FUN=function(x){unlist(strsplit(x,split=" +"))}))
		colnames(rules.mat) <- c("Parameter1_Name","relationship","Parameter2_Name")
	} else {
		rules.exist <- FALSE
		rules.mat <- NULL
	}
	if(length(est.parameters.lines)>0){
		parameters.exist  <- TRUE
		for(i in 1:length(variable.names)){
			est.parameters.lines <- gsub(variable.names[i],paste0(" ",variable.names[i]," "),est.parameters.lines,fixed=T)
		}
		est.parameters.lines <- gsub("^ +","",est.parameters.lines)
		est.parameters.lines <- as.list(gsub(" +$","",est.parameters.lines))
		parameters.mat    <- do.call(rbind,lapply(X=est.parameters.lines,FUN=function(x){unlist(strsplit(x,split=" +"))}))
		colnames(parameters.mat) <- c("IntegerThen1","ParameterName","Distribution","Min","Max","Output")
		### For each parameter, check that the "Max" value is greater than the "Min" value
		if(!all(as.numeric(parameters.mat[,"Max"]) > as.numeric(parameters.mat[,"Min"]))){
			offending.parameters <- paste(parameters.mat[!(as.numeric(parameters.mat[,"Max"]) > as.numeric(parameters.mat[,"Min"])),"ParameterName"],collapse=", ")
			stop.message <- paste0("In '.est' file: ",offending.parameters,": max value must be greater than min value")
			stop(stop.message)
		}
		### Check that names of distributions used for simple parameters are valid
		if(!all(parameters.mat[,"Distribution"] %in% c("unif","logunif"))){
			stop("simple parameters must have a prior distribution of type 'unif' or 'logunif'")
		}
		simple.pars  <- parameters.mat[,"ParameterName"]
	} else {
		parameters.exist <- FALSE
		parameters.mat   <- NULL
	}
	if(length(est.complex.lines)>0){
		if(!parameters.exist){
			stop("Simple parameters must be defined when complex parameters are defined")
		}
		# puts a space in front of and behind each parameter name
		for(i in 1:length(variable.names)){
			est.complex.lines <- gsub(variable.names[i],paste0(" ",variable.names[i]," "),est.complex.lines,fixed=T)
		}
		# puts a space in front of and behind each operator type. This will cause a problem if operators are used as part of parameter names.
		operator.types <- c("=","+","-","*","/","%min%","%max%","^","%abs%")
		for(i in 1:length(operator.types)){
			est.complex.lines <- gsub(operator.types[i],paste0(" ",operator.types[i]," "),est.complex.lines,fixed=T)
		}
		# removes spaces at beginning of strings
		est.complex.lines <- gsub("^ +","",est.complex.lines)
		# removes spaces at end of strings
		est.complex.lines <- as.list(gsub(" +$","",est.complex.lines))
		complex.exist <- TRUE
		# split by strings of spaces
		complex.list  <- lapply(X=est.complex.lines,FUN=function(x){unlist(strsplit(x,split=" +"))})
		# add empty entries to vectors that are shorter than the longest vector in complex.list, until all vectors are as long as the longest vector
		# The empty entries must be added before the last entry, which must be non-empty
		complex.list.lengths <- lengths(complex.list)
		if(any(complex.list.lengths    < max(complex.list.lengths))){
			complex.list.lengths.diffs <- max(complex.list.lengths) - complex.list.lengths
			for(i in 1:length(complex.list)){
				complex.list[[i]] <- c(complex.list[[i]], rep("",complex.list.lengths.diffs[i]))
			}
		}
		complex.mat   <- do.call(rbind,complex.list)
		if(ncol(complex.mat)>5){
			complex.mat.values <- complex.mat[,4:(ncol(complex.mat)-1),drop=F]
			complex.value <- do.call(rbind,lapply(X=1:nrow(complex.mat.values),FUN=function(x){paste(complex.mat.values[x,],collapse=" ")}))
			#complex.value <- as.matrix(paste(complex.mat[,4:(ncol(complex.mat)-1)],collapse=" "),byrow=T)
			complex.mat   <- cbind(complex.mat[,1:3,drop=F],complex.value,complex.mat[,ncol(complex.mat),drop=F])
		}
		colnames(complex.mat) <- c("IntegerThen1","ParameterName","function","Value","Output")
		complex.mat  <- gsub(" ","",complex.mat)
		complex.pars <- complex.mat[,"ParameterName"]
		all.pars     <- c(simple.pars,complex.pars)
		#### Making sure that parameters defined as complex are actually complex and not constants
		obj.mat <- matrix(rep(complex.mat[,"Value"],length(all.pars)),ncol=length(all.pars))
		pat.mat <- matrix(rep(all.pars,nrow(obj.mat)),ncol=length(all.pars),byrow=T)
		# Logical matrix indicating if the ith complex parameter uses the jth simple parameter.
		complex.of.all.mat           <- matrix(vapply(X=1:length(obj.mat),FUN=function(x){any(grep(pat.mat[x],obj.mat[x],fixed=T))},FUN.VALUE=TRUE),nrow=nrow(obj.mat))
		rownames(complex.of.all.mat) <- complex.pars
		colnames(complex.of.all.mat) <- all.pars
		complex.of.all               <- apply(X=complex.of.all.mat,MARGIN=1,FUN=any)
		if(!all(complex.of.all)){
			stop("The value of each complex parameter must be a function of at least one previously defined simple or complex parameter")
		}
	} else {
		complex.exist <- FALSE
		complex.mat   <- NULL
		complex.pars  <- NULL
	}
	initialization.attempt <- 1
	initialization.max     <- 100
	initialization.passed  <- FALSE
	while(initialization.attempt <= initialization.max & !initialization.passed){
		### Obtaining values to use for parameters.
		if(parameters.exist){
			parameters.meanval <- list(); length(parameters.meanval) <- nrow(parameters.mat)
			parameters.sampled <- list(); length(parameters.sampled) <- nrow(parameters.mat)
			for(i in 1:nrow(parameters.mat)){
				if(parameters.mat[i,"Distribution"]=="unif"){
					#### Keep the next line even though its commented out
					if(initialization.attempt==1 & use.mean){
						parameters.sampled[[i]] <- runif(100000,as.numeric(parameters.mat[i,"Min"]),as.numeric(parameters.mat[i,"Max"]))
						parameters.meanval[[i]] <- mean(c(as.numeric(parameters.mat[i,"Min"]),as.numeric(parameters.mat[i,"Max"])))
					} else {
						parameters.sampled[[i]] <- runif(1,as.numeric(parameters.mat[i,"Min"]),as.numeric(parameters.mat[i,"Max"]))
						parameters.meanval[[i]] <- parameters.sampled[[i]]
					}
				}
				if(parameters.mat[i,"Distribution"]=="logunif"){
					if(initialization.attempt==1 & use.mean){
						parameters.sampled[[i]] <- rlogunif(100000,as.numeric(parameters.mat[i,"Min"]),as.numeric(parameters.mat[i,"Max"]))
						parameters.meanval[[i]] <- mean(parameters.sampled[[i]])
					} else {
						parameters.sampled[[i]] <- rlogunif(1,as.numeric(parameters.mat[i,"Min"]),as.numeric(parameters.mat[i,"Max"]))
						parameters.meanval[[i]] <- parameters.sampled[[i]]
					}
				}
				if(parameters.mat[i,"IntegerThen1"]==1){
					parameters.sampled[[i]] <- round(parameters.sampled[[i]])
					parameters.meanval[[i]] <- round(parameters.meanval[[i]])
				}
			}
			parameters.meanval        <- unlist(parameters.meanval)
			names(parameters.meanval) <- parameters.mat[,"ParameterName"]
			parameters.sampled.mat    <- do.call(cbind,parameters.sampled)
			colnames(parameters.sampled.mat) <- parameters.mat[,"ParameterName"]
			### Holding each distributions as an object, using the 'distr' package
		}
		
		if(complex.exist){
			for(i in 1:length(all.pars)){
				if(i == 1){
					complex.operators   <- gsub(all.pars[i],"",complex.mat[,"Value"],fixed=T)
				} else {
					complex.operators   <- gsub(all.pars[i],"",complex.operators,fixed=T)
				}
			}
			for(i in 0:9){
				complex.operators   <- gsub(i,"",complex.operators,fixed=T)
			}
			names(complex.operators) <- complex.mat[,"ParameterName"]
			#print(complex.operators)
			#### This may bring a problem for abs or exp parameters
			complex.operands <- list(); length(complex.operands) <- length(complex.operators)
			for(i in 1:length(complex.operators)){
				#if(i == 1){
					if(complex.operators[i] %in% c("log()","abs()","exp()")){
						complex.operands[[i]] <- c("","")
						next
					}
					
					complex.operands[[i]]   <- gsub(" ","",unlist(strsplit(complex.mat[i,"Value"],split=complex.operators[i], fixed=T)))
					#complex.operands[[i]]   <- unlist(sapply(1:length(complex.mat[,"Value"]),FUN=function(x){unlist(strsplit(complex.mat[x,"Value"],split=complex.operators[i],fixed=T))}))
					#complex.rightoperand <- sapply(1:length(complex.mat[,"Value"]),FUN=function(x){unlist(strsplit(complex.mat[x,"Value"],split=complex.operators[i],fixed=T))})
					#complex.leftoperand  <- sapply(1:length(complex.mat[,"Value"]),FUN=function(x){unlist(strsplit(complex.mat[x,"Value"],split=complex.operators[i],fixed=T))[1]})
					#complex.rightoperand <- sapply(1:length(complex.mat[,"Value"]),FUN=function(x){unlist(strsplit(complex.mat[x,"Value"],split=complex.operators[i],fixed=T))[2]})
				#} else {
				#	complex.operands[[i]]   <- unlist(sapply(1:length(complex.operands),FUN=function(x){unlist(strsplit(complex.operands[x],split=complex.operators[i],fixed=T))}))
					#complex.rightoperand <- sapply(1:length(complex.rightoperand),FUN=function(x){unlist(strsplit(complex.rightoperand[x],split=complex.operators[i],fixed=T))})
				#}
			}
			complex.operands.mat            <- do.call(rbind,complex.operands)
			#complex.operands               <- gsub(" ","",complex.operands)
			#complex.operands.mat           <- matrix(complex.operands,ncol=2,byrow=TRUE)
			if(!is.null(complex.operands.mat)){
				colnames(complex.operands.mat) <- c("left","right")
				rownames(complex.operands.mat) <- names(complex.operators)[lengths(complex.operands)>0]
				#print(complex.operands.mat)
				#test.left  <- complex.operands.mat[,"left"]
				#test.right <- complex.operands.mat[,"right"]
				complex.leftoperand  <- complex.operands.mat[,"left"]
				complex.rightoperand <- complex.operands.mat[,"right"]
				complex.value.strings <- paste(complex.leftoperand,complex.operators[lengths(complex.operands)>0],complex.rightoperand,sep=" ")
			} else {
				complex.leftoperand  <- NULL
				complex.rightoperand <- NULL
			}
			#print(complex.value.strings)
			#print(complex.mat)
			#print(test.right)
			#print(complex.operands.mat[,"left"])
			#print(complex.operands.mat[,"right"])
			complex.val.string  <- list(); length(complex.val.string)  <- nrow(complex.mat)
			complex.text.string <- list(); length(complex.text.string) <- nrow(complex.mat)
			complex.val <- list(); length(complex.val) <- nrow(complex.mat)
			names(complex.val)         <- complex.mat[,"ParameterName"]
			names(complex.val.string)  <- complex.mat[,"ParameterName"]
			names(complex.text.string) <- complex.mat[,"ParameterName"]
			for(i in 1:nrow(complex.mat)){
				if(complex.operators[i] %in% c("log()","abs()","exp()")){
					value.string      <- complex.mat[i,"Value"]
				} else {
					value.string      <- complex.value.strings[i]
				}
				value.vector      <- unlist(strsplit(value.string,split=" "))
				search.parameters <- unlist(lapply(X=1:nrow(parameters.mat),FUN=function(z){any(grep(names(parameters.meanval)[z],value.vector,fixed=T))}))
				if(any(search.parameters)){
					parameters.temp         <- names(parameters.meanval)[which(search.parameters)]
					parameters.meanval.temp <- parameters.meanval[which(search.parameters)]
				} else {
					parameters.temp         <- NULL
					parameters.meanval.temp <- NULL
				}
				if(i>1){
					### look in previously evaluated complex parameters
					search.complex <- unlist(lapply(X=1:(i-1),FUN=function(z){any(grep(names(complex.val)[z],value.vector,fixed=T))}))
					if(any(search.complex)){
						complex.temp     <- names(complex.val)[which(search.complex)]
						complex.val.temp <- complex.val[which(search.complex)]
					} else {
						complex.temp     <- NULL
						complex.val.temp <- NULL
					}
				} else {
					complex.temp     <- NULL
					complex.val.temp <- NULL
				}
				value.vector2       <- value.vector
				value.string2       <- value.string
				text.string         <- value.string
				parameters.temp     <- c(parameters.temp,complex.temp)
				parameters.val.temp <- c(parameters.meanval.temp,complex.val.temp)
				if(!is.null(parameters.temp)){
					for(j in 1:length(parameters.temp)){
						value.string2 <- gsub(parameters.temp[j],parameters.val.temp[j],value.string2,fixed=T)
					}
				}
				if(any(grep("%min%",value.string2,fixed=T))){
					complex.val[[i]]         <- min(unlist(parameters.val.temp))
					complex.val.string[[i]]  <- paste0("min(",parameters.val.temp[1],",",parameters.val.temp[2],")")
					complex.text.string[[i]] <- paste0("min(",names(parameters.val.temp)[1],",",names(parameters.val.temp)[2],")")
				} else {
					if(any(grep("%max%",value.string2,fixed=T))){
						complex.val[[i]]         <- max(unlist(parameters.val.temp))
						complex.val.string[[i]]  <- paste0("max(",parameters.val.temp[1],",",parameters.val.temp[2],")")
						complex.text.string[[i]] <- paste0("max(",names(parameters.val.temp)[1],",",names(parameters.val.temp)[2],")")
					} else {
						complex.val.string[[i]]  <- value.string2
						complex.val[[i]]         <- eval(parse(text=value.string2))
						complex.text.string[[i]] <- text.string
					}
				}
			}
			complex.val                <- unlist(complex.val)
			names(complex.val)         <- complex.mat[,"ParameterName"]
			complex.val.string         <- unlist(complex.val.string)
			names(complex.val.string)  <- complex.mat[,"ParameterName"]
			complex.text.string        <- unlist(complex.text.string)
			names(complex.text.string) <- complex.mat[,"ParameterName"]
		}
		### If simple and complex parameters exist
		if(parameters.exist & complex.exist){
			variables <- c(parameters.meanval,complex.val)
		}
		### If all parameters are simple
		if(parameters.exist & !complex.exist){
			variables <- c(parameters.meanval)
		}
		#### Check that the rules are obeyed. If not obeyed, resample from prior until rules are satisfied, or until maximum initialization attempts reached.
		if(rules.exist){
			variables.with.rules <- variables[names(variables) %in% rules.mat]
			# initialize rules.values.mat by setting it equal to rules.mat
			rules.values.mat <- rules.mat
			for(i in 1:length(variables.with.rules)){
				rules.values.mat[rules.values.mat == names(variables.with.rules)[i]] <- variables.with.rules[i]
			}
			rules.vals.strings <- apply(X=rules.values.mat,MARGIN=1,FUN=paste,collapse=" ")
			#print(rules.vals.strings)
			rules.evaluated    <- sapply(X=1:length(rules.vals.strings), FUN=function(x){eval(parse(text=rules.vals.strings[x]))})
			
			if(!all(rules.evaluated)){
				if(initialization.attempt  < initialization.max){
					initialization.attempt <- (initialization.attempt + 1)
					#initialization.passed  <- FALSE
					initialization.check1 <- FALSE
					next
				} else {
					stop(paste("Rules",which(!rules.evaluated),"failed to initialize as true. Model may be valid but unable to find solution for plotting."))
				}
			} else {
				#initialization.attempt <- initialization.max+1
				initialization.check1  <- TRUE
			}
		} else {
			initialization.check1 <- TRUE
		}
		### vector with lines containing historical events
		events.lines.temp    <- tpl.lines[(grep(pattern="^//historical event:", tpl.lines)+2):(grep(pattern="^//Number of independent loci", tpl.lines)-1)]
		### Checks if end-of-line comments exist, and warns if they do because fsc26 seems to fail when c++ comments are in unusual places
		if(any(grep("\\\\",events.lines.temp)) | any(grep("\\/\\/",events.lines.temp)) | any(grep("#",events.lines.temp))){
			events.lines.temp <- gsub("\\\\.+","",events.lines.temp)
			events.lines.temp <- gsub("\\/\\/.+","",events.lines.temp)
			events.lines.temp <- gsub("#.+","",events.lines.temp)
			if(warn.about.comments){
				warning("Comments were found on one or more lines defining historical events in '.tpl' file; these comments may not be supported by fsc.")
			}
		}
		### Remove whitespace before linebreaks
		events.lines.temp    <- gsub(" +$","",events.lines.temp)
		### Hold events as a list
		events.lines         <- as.list(events.lines.temp)
		num.events           <- as.numeric(unlist(strsplit(tpl.lines[(grep(pattern="^//historical event:", tpl.lines)+1)],split=" "))[1])
		events.mat           <- do.call(rbind,lapply(X=events.lines,FUN=function(x){unlist(strsplit(x,split=" "))}))
		colnames(events.mat) <- c("time","source","sink","migrants","newsize","growthrate","migration matrix")
		if(any(names(variables) %in% events.mat)){
			variables.in.events  <- variables[names(variables) %in% events.mat]
			### version of events.mat that contains values instead of variable names
			events.mat.values <- events.mat
			for(i in 1:length(variables.in.events)){
				variable.temp     <- variables.in.events[i]
				events.mat.values <- gsub(names(variable.temp),variable.temp,events.mat.values,fixed=T)
			}
		} else {
			events.mat.values <- events.mat
		}
		# Sort rows of events.mat.values and events.mat by increasing value in the first column and then by presence or absence of "keep" for any value in a row
		rows.reorder0     <- rep(0,nrow(events.mat.values))
		rows.reorder0[sapply(X=1:nrow(events.mat.values),FUN=function(x){"keep" %in% events.mat.values[x,]})] <- 1
		rows.reorder      <- order( as.numeric(events.mat.values[,"time"]), rows.reorder0)
		events.mat.values <- events.mat.values[rows.reorder,,drop=F]
		events.mat        <- events.mat[rows.reorder,,drop=F]
		### Number of populations at present (i.e., zero generations ago) and their population sizes
		npops.present    <- as.numeric(tpl.lines[2])
		
		### names of variables that specify intial population sizes
		popsizes.initial <- tpl.lines[(grep("//Population effective sizes",tpl.lines)+1):(grep("//Population effective sizes",tpl.lines)+npops.present)]
		if(any(grep("^[0-9]+$",popsizes.initial))){
			popsizes.pars    <- popsizes.initial[!grep("^[0-9]+$",popsizes.initial)]
		} else {
			popsizes.pars    <- popsizes.initial
		}
		if(length(popsizes.pars)>0){
			popsizes.present <- parameters.meanval[names(parameters.meanval) %in%  popsizes.pars]
		} else {
			popsizes.present <- as.numeric(popsizes.initial)
		}
		### names of variables that specify changes in growth rate
		if(any(is.na(suppressWarnings(as.numeric(events.mat[,"growthrate"]))))){
			growthrate.pars <- events.mat[is.na(suppressWarnings(as.numeric(events.mat[,"growthrate"]))),"growthrate"]
		} else {
			growthrate.pars <- NULL
		}
		### names of variables that specify instantaneous changes in size
		if(any(is.na(suppressWarnings(as.numeric(events.mat[,"newsize"]))))){
			newsize.pars <- events.mat[is.na(suppressWarnings(as.numeric(events.mat[,"newsize"]))),"newsize"]
		} else {
			newsize.pars <- NULL
		}
		### names of variables that specify instantaneous migration
		if(any(is.na(suppressWarnings(as.numeric(events.mat[,"migrants"]))))){
			migrants.pars <- events.mat[is.na(suppressWarnings(as.numeric(events.mat[,"migrants"]))),"migrants"]
		} else {
			migrants.pars <- NULL
		}
		### Number of migration matrices
		num.mig.matrices <- as.numeric(tpl.lines[grep("//Number of migration matrices",tpl.lines)+1])
		# Replace "keep" values with their migration matrix index
		if(any(events.mat.values[,"migration matrix"]=="keep")){
			for(i in 1:nrow(events.mat.values)){
				if(i==1 & events.mat.values[i,"migration matrix"]=="keep"){
					events.mat.values[i,"migration matrix"] <- "0"
				}
				if(i!=1 & events.mat.values[i,"migration matrix"]=="keep"){
					events.mat.values[i,"migration matrix"] <- events.mat.values[(i-1),"migration matrix"]
				}
			}
		}
		## Initial values for growth rate
		initial.growthrate <- as.numeric(tpl.lines[(grep("//Growth rates",tpl.lines)+1):(grep("//Growth rates",tpl.lines)+npops.present)])
		### If 'keep' is used as a directive for growth rate, remove the growthrates column from the historical events matrix so that this information can be stored in a list of numerical vectors
		growthrates.list <- list(); length(growthrates.list) <- (nrow(events.mat.values)+1)
		growthrates.list[[1]] <- as.numeric(initial.growthrate)
		for(i in 1:nrow(events.mat.values)){
			if(events.mat.values[i,"growthrate"]=="keep"){
				growthrates.list[[i+1]] <- growthrates.list[[i]]
			} else {
				### Not sure how to incorporate popsize change information yet
				growthrates.list[[i+1]] <- as.numeric(rep(events.mat.values[i,"growthrate"],npops.present))
			}
		}
		# Remove growthrates column from events.mat.values now that this info has been moved to a list.
		events.mat.values <- events.mat.values[, !(colnames(events.mat.values) %in% "growthrate"),drop=F]
		### change the mode of events.mat.values to numeric
		mode(events.mat.values) <- "numeric"
		time.list <- list(); length(time.list) <- nrow(events.mat.values)+1
		for(i in 1:length(time.list)){
			if(i == 1){
				start.time <- 0
			} else {
				start.time <- events.mat.values[i-1,"time"]
			}
			if(i < length(time.list)){
				end.time   <- events.mat.values[i,"time"]-1
			} else {
				end.time   <- start.time*1.5
			}
			time.list[[i]] <- rbind(start.time,end.time)
		}
		### Number of time periods
		num.periods <- nrow(events.mat.values)+1
		#}
		### Internal names to use for populations
		popnames.short <- paste0("pop",c(0:(length(popsizes.present)-1)))
		### Don't delete period.sizes.v2
		period.sizes.v2 <- list(); length(period.sizes.v2) <- num.periods
		for(i in 1:length(period.sizes.v2)){
			if(i == 1){
				start.sizes <- popsizes.present
				names(start.sizes) <- popnames.short
				GROWTH      <- start.sizes*growthrates.list[[i]]*(time.list[[i]]["end.time",]-time.list[[i]]["start.time",])
				end.sizes   <- start.sizes + GROWTH
				period.sizes.v2[[i]] <- rbind(start.sizes,end.sizes)
			} else {
				sourcepop   <- events.mat.values[i-1,"source"]
				sinkpop     <- events.mat.values[i-1,"sink"]
				MIGRANTS    <- period.sizes.v2[[i-1]]["end.sizes",sourcepop+1] * events.mat.values[i-1,"migrants"]
				# temporarily set start sizes equal to previous end.sizes
				start.sizes <- period.sizes.v2[[i-1]]["end.sizes",]
				start.sizes[sourcepop+1] <- start.sizes[sourcepop+1] - MIGRANTS
				start.sizes[sinkpop+1]   <- start.sizes[sinkpop+1] + MIGRANTS
				# Now modify by the value of newsize
				start.sizes <- start.sizes * events.mat.values[i-1,"newsize"]
				GROWTH      <- start.sizes * growthrates.list[[i]]*(time.list[[i]]["end.time",]-time.list[[i]]["start.time",])
				end.sizes   <- start.sizes + GROWTH
				period.sizes.v2[[i]] <- rbind(start.sizes,end.sizes)
			}
		}
		period.sizes.v3 <- list(); length(period.sizes.v3) <- length(period.sizes.v2)
		for(i in 1:length(period.sizes.v3)){
			period.sizes.v3[[i]] <- cbind(time.list[[i]],period.sizes.v2[[i]])
		}
		## list with migration matrices
		if(num.mig.matrices>0){
			matrix.header.lines <- grep("//migrationmatrix", tpl.lines)
			migration.matrices  <- lapply(X=1:length(matrix.header.lines),FUN=function(x){matrix(unlist(strsplit(tpl.lines[(matrix.header.lines[x]+1):(matrix.header.lines[x]+npops.present)],split=" ")),ncol=npops.present,nrow=npops.present,byrow=T)})
			for(i in 1:length(migration.matrices)){
				colnames(migration.matrices[[i]]) <- popnames.short
				rownames(migration.matrices[[i]]) <- popnames.short
			}
		} else {
			migration.matrices  <- NULL
		}
		slices.popsizes.all <- do.call(rbind,period.sizes.v3)[,-1]
		### Unique combinations of presence/absence among populations at any point in time.
		pops.present <- slices.popsizes.all
		pops.present[pops.present > 0] <- 1
		unique.pops.present <- unique(pops.present)
		if(any(unique.pops.present < 0)){
			if(initialization.attempt < initialization.max){
				initialization.attempt <- initialization.attempt + 1
				next
			} else {
				stop("Model may be valid, but plotting failed when using means of priors because negative population sizes produced. A future version of this function may resample prior to find valid solution.")
			}
		} else {
			initialization.check2 <- TRUE
		}
		if(initialization.check1 & initialization.check2){
			initialization.passed  <- TRUE
			if(initialization.attempt>1){
				print(paste("initialized after",initialization.attempt,"attempts"))
			}
		}
	}

	max.popsize.eachpop.list <- list(); length(max.popsize.eachpop.list) <- nrow(unique.pops.present)
	names(max.popsize.eachpop.list) <- apply(X= unique.pops.present,MARGIN=1,FUN=paste,collapse=" ")
	for(i in 1:nrow(unique.pops.present)){
		combo.temp <- unique.pops.present[i,]
		## subset of slices.popsizes.all containing columns with the ith combination of pops and rows without zeros.
		slices.popsizes.i <- slices.popsizes.all[,which(unique.pops.present[i,]==1)]
		rows.keep <- apply(X=pops.present,MARGIN=1,FUN=function(x){all(x==combo.temp)})
		slices.popsizes.temp <- slices.popsizes.all[rows.keep,]
		max.popsize.eachpop.temp <- apply(X=slices.popsizes.temp,MARGIN=2,FUN=max)
		max.popsize.eachpop.list[[i]] <- max.popsize.eachpop.temp[max.popsize.eachpop.temp>0]
	}
	####
	y.extents <- list(); length(y.extents) <- length(max.popsize.eachpop.list)
	names(y.extents) <- names(max.popsize.eachpop.list)
	for(z in 1:length(max.popsize.eachpop.list)){
		ytops <- list(); length(ytops) <- length(max.popsize.eachpop.list[[z]])
		for(i in 1:length(ytops)){
			if(i==1){
				top.current <- 0
			}
			ytops[[i]] <- max.popsize.eachpop.list[[z]][i] + (top.current*pops.spacer)
			top.current <- ytops[[i]]
		}
		ytops    <- unlist(ytops)
		ybottoms <- ytops-max.popsize.eachpop.list[[z]]
		ymids    <- colMeans(rbind(ybottoms,ytops))
		y.extents[[z]] <- rbind(ybottoms,ymids,ytops)
	}
	shapes.list <- list(); length(shapes.list) <- length(period.sizes.v3)
	for(i in 1:length(period.sizes.v3)){
		# matrix with popsizes and start and end of ith time period
		period.mat     <- period.sizes.v3[[i]]
		# number of populations in the ith period
		period.numpops      <- length(which(apply(X=period.mat[,-1],MARGIN=2,FUN=function(x){all(x!=0)})))
		# number of populations in each of the first through ith periods
		if(i==1){
			period.numpops.vector <- period.numpops
		} else {
			period.numpops.vector <- c(period.numpops.vector,period.numpops)
		}
		### y-axis midlines for populations
		pops.present.temp <- paste(as.numeric(unique(period.mat[,-1] > 0)),collapse=" ")
		y.extents.temp    <- y.extents[[which(names(y.extents) %in% pops.present.temp)]]
		period.mat        <- period.mat[,c(1,which(colnames(period.mat) %in% colnames(y.extents.temp))),drop=F]
		shapes.list.i     <- list(); length(shapes.list.i) <- ncol(y.extents.temp)
		for(j in 1:ncol(y.extents.temp)){
			## start and end population size for population j-1 during period i
			period.mat.j       <- period.mat[,c(1,(j+1))]
			x_left             <- period.mat.j[1,1]
			x_right            <- period.mat.j[2,1]
			y_mid              <- y.extents.temp["ymids",j]
			y_top_left         <- y_mid+(period.mat.j[1,2]/2)
			y_top_right        <- y_mid+(period.mat.j[2,2]/2)
			y_bottom_left      <- y_mid-(period.mat.j[1,2]/2)
			y_bottom_right     <- y_mid-(period.mat.j[2,2]/2)
			shape.points.temp  <- c(period=i,x_left,x_right,y_top_left,y_top_right,y_bottom_left,y_bottom_right)
			shapes.list.i[[j]] <- shape.points.temp
		}
		shapes.list[[i]] <- do.call(rbind,shapes.list.i)
		colnames(shapes.list[[i]]) <- c("period","x_left","x_right","y_top_left","y_top_right","y_bottom_left","y_bottom_right")
		rownames(shapes.list[[i]]) <- colnames(period.mat)[-1]
	}
	if(npops.present>1){
		# Check if at least one divergence/coalescent event occurs, and if so create matrices holding info on coalescence and divergence events, and possibly also fusions
		if(any(events.mat.values[,"migrants",drop=F]==1)){
			### Rows of historical events matrix that define divergences (coalescences)
			coalescence.event.rows <- which(events.mat.values[,"migrants",drop=F]==1)
			### Matrix with coalescence events
			coalescence.events.mat <- events.mat.values[coalescence.event.rows,,drop=F]
			### reverse-order of divergence matrix
			divergence.events.mat <- coalescence.events.mat[rev(1:nrow(coalescence.events.mat)),,drop=F]
			### matrix with only the source and sink columns of the divergence.events.mat
			divergence.mat <- divergence.events.mat[,c("source","sink"),drop=F]
			### populations that were sources during at least one coalescence event
			coal.sources.unique <- unique(divergence.mat[,"source"])
			### Historical events in which a nonzero fraction of a population is transferred from between two different populations.
			migration.events.rows <- which(events.mat.values[,"migrants"]!=0 & events.mat.values[,"source"] != events.mat.values[,"sink"])
			migration.events.mat  <- events.mat.values[migration.events.rows,,drop=F]
			## the length of fusion.events mat is initially set to the number of divergence events
			fusion.events <- list(); length(fusion.events) <- nrow(divergence.mat)
			for(i in 1:length(coal.sources.unique)){
				# divergence events in which the ith population of coal.sources.unique is the source population
				divergence.mat.i <- divergence.events.mat[which(divergence.events.mat[,"source"] == coal.sources.unique),,drop=F]
				# a list to hold only the fusions reversing divergences with the ith population as a source
				fusion.events.i <- list(); length(fusion.events.i) <- nrow(divergence.mat.i)
				for(j in 1:nrow(divergence.mat.i)){
					divergence.ij.time <- divergence.mat.i[j,"time"]
					## check if any more recent migration events involve the ith source pop as a sink pop
					if(any(migration.events.mat[,"sink"]==coal.sources.unique[i] & migration.events.mat[,"time"] < divergence.ij.time)){
						nextevent.row <- which(migration.events.mat[,"sink"]==coal.sources.unique[i] & migration.events.mat[,"time"] < divergence.ij.time)
						## the next event should be a fusion event
						nextevent     <- migration.events.mat[nextevent.row,,drop=F]
						fusion.events.i[[j]] <- nextevent
					} else {
						fusion.events.i[[j]] <- rep(NA,6)
					}
				}
				fusion.events.i <- do.call(rbind,fusion.events.i)
				fusion.events[[i]] <- fusion.events.i
			}
			fusion.events.mat <- do.call(rbind,fusion.events)
			if(all(is.na(fusion.events.mat[,1]))){
				fusions.exist=F
			} else {
				fusions.exist=T
				### Remove NA rows that were created for divergences that not followed by fusions
				if(any(is.na(fusion.events.mat[,1]))){
						fusion.events.mat <- fusion.events.mat[-which(is.na(fusion.events.mat[,1])),,drop=F]
				}
				### Ensure that fusion events are ordered from most ancient to most recent.
				fusion.events.mat[order(fusion.events.mat[,"time"],decreasing=T),,drop=F]
				### Simpler version of fusion.events.mat containing only the source and sink columns. This will be used during plotting.
				fusions.mat <- fusion.events.mat[,c("source","sink"),drop=F]
			}
		}
		### Construct the newick string representing relationships among populations, and then convert to phylo. Used to test if priors are consistent with a bifurcating phylogeny.
		for(i in 1:nrow(divergence.events.mat)){
			divergence.i           <- divergence.events.mat[i,c("source","sink"),drop=F]
			divergence.string.temp <- paste0("(",paste0(divergence.i,collapse=","),")")
			
			if(i==1){
				divergence.string  <- divergence.string.temp
			} else {
				sink.i <- divergence.i[,"sink"]
				### check if the current sink is in the divergence string and if so replace the sink value in the divergence string with the value of divergence.string.temp
				### This works as long as you dont have more than 10 populations (pop0-pop9)
				if(any(grep(sink.i,divergence.string))){
					divergence.string <- gsub(sink.i,divergence.string.temp,divergence.string)
				}
			}
			if(i==nrow(divergence.events.mat)){
				divergence.string <- paste0("(",divergence.string,");")
			}
		}
	#	tree.phylo <- ape::read.tree(text=divergence.string)
	#	plot(tree.phylo)
	#	result.test <- recordPlot()
	#	return(result.test)
		if(length(shapes.list)>1){
			# initialize object to adjust during loop
			shapes.list.v2 <- rev(shapes.list)
			for(i in 2:length(shapes.list.v2)){
				# initial settings for counters. Will add +1 to div.temp after each divergence event, and +to fus.temp after each fusion event
				if(i==2){
					div.temp=0
					fus.temp=0
				}
				shapes.list.v2.temp <- shapes.list.v2[[i]]
				### Checks if all populations in the current period were also in the previous period. If so, it means no divergence events in that period.
				if(length(setdiff(rownames(shapes.list.v2.temp),rownames(shapes.list.v2[[i-1]])))==0){
					### Checks if all populations in the previous period are also in the current period. If so, it means no fusion events in that period.
					if(length(setdiff(rownames(shapes.list.v2[[i-1]]), rownames(shapes.list.v2.temp)))==0){
						# Align all popshapes with their previous period popshapes
						pops.other    <- rownames(shapes.list.v2.temp)
						oldmidlines.i <- apply(shapes.list.v2.temp[,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right"),drop=F],MARGIN=1,FUN=mean)
						newmidlines.i <- apply(shapes.list.v2[[i-1]][,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right"),drop=F],MARGIN=1,FUN=mean)
						dif.midlines  <- newmidlines.i-oldmidlines.i
						### updates shapes.list.v2
						pops.other.vals <- shapes.list.v2[[i]][,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right"),drop=F] + dif.midlines
						shapes.list.v2[[i]][,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right")] <- pops.other.vals
					} else {
						fus.temp = fus.temp + 1
						pop.lost <- setdiff(rownames(shapes.list.v2[[i-1]]), rownames(shapes.list.v2.temp))
						# Population that was sister to pop.lost in the previous period.
						pop.sis  <- paste0("pop",fusions.mat[fus.temp,!(fusions.mat[fus.temp,] %in% gsub("pop","",pop.lost))])
						### current (not adjusted) midline of pop.sis in the current period
						oldmidlines.pop.sis  <- apply(shapes.list.v2.temp[pop.sis,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right"),drop=F],MARGIN=1,FUN=mean)
						### Get the mean y value of y_bottom_left of the upper shape and y_top_left of the lower shape for pop.lost and pop.sis in the previous generation
						if(shapes.list.v2[[i-1]][pop.sis,"y_bottom_left",drop=F] > shapes.list.v2[[i-1]][pop.lost,"y_top_left",drop=F]){
							newmidlines.pop.sis <- mean(c(shapes.list.v2[[i-1]][pop.sis,"y_bottom_left",drop=F],shapes.list.v2[[i-1]][pop.lost,"y_top_left",drop=F]))
						} else {
							newmidlines.pop.sis <- mean(c(shapes.list.v2[[i-1]][pop.lost,"y_bottom_left",drop=F],shapes.list.v2[[i-1]][pop.sis,"y_top_left",drop=F]))
						}
						### Amount and direction to shift y values of pop.sis
						diff.sis <- newmidlines.pop.sis - oldmidlines.pop.sis
						### New y values for pop.sis
						pop.sis.vals <- shapes.list.v2[[i]][pop.sis,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right"),drop=F] + c(diff.sis)
						if(length(rownames(shapes.list.v2.temp))>1){
							pops.other    <- setdiff(rownames(shapes.list.v2.temp),c(pop.sis))
							oldmidlines.i <- apply(shapes.list.v2.temp[pops.other,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right"),drop=F],MARGIN=1,FUN=mean)
							newmidlines.i <- apply(shapes.list.v2[[i-1]][pops.other,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right"),drop=F],MARGIN=1,FUN=mean)
							diff.midlines <- newmidlines.i-oldmidlines.i
							### updates shapes.list.v2 midlines for shapes not involved in a divergence event
							pops.other.vals <- shapes.list.v2[[i]][pops.other,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right"),drop=F] + diff.midlines
							shapes.list.v2[[i]][pops.other,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right")] <- pops.other.vals
						}
					}
				} else {
					div.temp = div.temp + 1
					# Determine which population is 'new', and call this pop.new
					pop.new <- setdiff(rownames(shapes.list.v2.temp),rownames(shapes.list.v2[[i-1]]))
					# Consult the divergence matrix to determine the sister population to pop.new
					pop.sis <- paste0("pop",divergence.mat[div.temp,!(divergence.mat[div.temp,] %in% gsub("pop","",pop.new))])
					### Want to set the mean of the daughter shape midlines equal to the parent shape midline.
					oldmidlines.pop.new <- apply(shapes.list.v2.temp[pop.new,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right"),drop=F],MARGIN=1,FUN=mean)
					oldmidlines.pop.sis <- apply(shapes.list.v2.temp[pop.sis,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right"),drop=F],MARGIN=1,FUN=mean)
					### Setting the midlines for the upper and lower shapes of period i equal to the y values of the top left and bottom left corners of the parent shape, respectively
					if(oldmidlines.pop.new > oldmidlines.pop.sis){
						newmidlines.pop.new <- shapes.list.v2[[i-1]][pop.sis,"y_top_left",drop=F]
						newmidlines.pop.sis <- shapes.list.v2[[i-1]][pop.sis,"y_bottom_left",drop=F]
					} else {
						newmidlines.pop.sis <- shapes.list.v2[[i-1]][pop.sis,"y_top_left",drop=F]
						newmidlines.pop.new <- shapes.list.v2[[i-1]][pop.sis,"y_bottom_left",drop=F]
					}
					diff.new <- newmidlines.pop.new - oldmidlines.pop.new
					diff.sis <- newmidlines.pop.sis - oldmidlines.pop.sis
			#		### updates shapes.list.v2 midlines for the two shapes that are involved in a divergence event
					pop.new.vals <- shapes.list.v2[[i]][pop.new,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right"),drop=F] + c(diff.new)
					pop.sis.vals <- shapes.list.v2[[i]][pop.sis,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right"),drop=F] + c(diff.sis)
					shapes.list.v2[[i]][pop.new,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right")] <- pop.new.vals
					shapes.list.v2[[i]][pop.sis,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right")] <- pop.sis.vals
					if(length(rownames(shapes.list.v2.temp))>2){
						pops.other    <- setdiff(rownames(shapes.list.v2.temp),c(pop.new,pop.sis))
						oldmidlines.i <- apply(shapes.list.v2.temp[pops.other,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right"),drop=F],MARGIN=1,FUN=mean)
						newmidlines.i <- apply(shapes.list.v2[[i-1]][pops.other,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right"),drop=F],MARGIN=1,FUN=mean)
						diff.midlines <- newmidlines.i-oldmidlines.i
						### updates shapes.list.v2 midlines for shapes not involved in a divergence event
						pops.other.vals <- shapes.list.v2[[i]][pops.other,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right"),drop=F] + diff.midlines
						shapes.list.v2[[i]][pops.other,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right")] <- pops.other.vals
					}
				}
			}
		}
		shapes.list <- shapes.list.v2
	}
	### adjust all y-coordinates to ensure that all parts of the plot are visible.
	shapes.mat <- do.call(rbind,shapes.list)
	if(min(shapes.mat[,c("y_bottom_left","y_bottom_right")]) < 0){
		y.adjust <- abs(min(shapes.mat[,c("y_bottom_left","y_bottom_right")]))
		shapes.mat[,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right")] <- shapes.mat[,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right")] + y.adjust
	}
	### Subsetting parameters names by type. Maybe move this to a slightly earlier spot.
	if(length(intersect(events.mat[,"time"],names(variables)))!=0){
		time.pars  <- intersect(events.mat[,"time"],names(variables))
	} else {
		time.pars <- NULL
	}
	#if(length(intersect(events.mat[,"time"],names(variables)))!=0){
	#	time.pars  <- intersect(events.mat[,"time"],names(variables))
	#}
	if(length(popsizes.pars)==0){
		popsizes.pars <- NULL
	}
	if(length(migration.matrices)!=0){
		migrationrate.pars <- intersect(unique(unlist(migration.matrices)),names(variables))
		if(length(migrationrate.pars)==0){
			migrationrate.pars <- NULL
		}
	} else {
		migrationrate.pars <- NULL
	}
	other.pars         <- setdiff(names(variables),c(time.pars,migrationrate.pars,popsizes.pars,migrants.pars,newsize.pars,growthrate.pars))
	if(length(other.pars)==0){
		other.pars <- NULL
		other.simple.pars  <- NULL
	} else {
		other.simple.pars  <- intersect(other.pars,simple.pars)
	}
	if(length(other.simple.pars)==0){
		other.simple.pars  <- NULL
	}
	### Separating variables.sampled.mat by parameters that are uniformly distributed and those that are loguniformly distributed
	# Not sure what to do if a complex parameter is evaluated from both a uniform and loguniformly distributed parameter
#	unif.simple.pars    <- intersect(simple.pars,parameters.mat[parameters.mat[,"Distribution"] %in% "unif","ParameterName"])
#	logunif.simple.pars <- intersect(simple.pars,parameters.mat[parameters.mat[,"Distribution"] %in% "logunif","ParameterName"])
	
	# The purpose of this is to classify the complex parameters as having either a uniform or loguniform probability density.
	if(length(complex.pars)>0){
			if(!is.null(complex.leftoperand) & !is.null(complex.rightoperand)){
			complex.dist.list <- list(); length(complex.dist.list) <- nrow(complex.mat)
			names(complex.dist.list) <- complex.mat[,"ParameterName"]
			for(i in 1:nrow(complex.mat)){
				if(complex.leftoperand[i] %in% simple.pars & any(grep("^[0-9]+$",complex.rightoperand[i]))){
					if(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Distribution"] == "logunif"){
						complex.dist.list <- NULL
					} else {
						if(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Distribution"] == "unif"){
							### uniform simple parameter divided by a constant
							if(complex.operators[i] == "/"){
								complex.min.temp       <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Min"])/as.numeric(complex.rightoperand[i])
								complex.max.temp       <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Max"])/as.numeric(complex.rightoperand[i])
								result.i               <- matrix(c(complex.mat[i,"ParameterName"],complex.min.temp,complex.max.temp,"Unif/Constant"),nrow=1,ncol=4)
								colnames(result.i)     <- c("Parameter","Min","Max","DistClass")
								complex.dist.list[[i]] <- as.data.frame(result.i)
							}
							### uniform simple parameter multiplied by a constant
							if(complex.operators[i] == "*"){
								complex.min.temp       <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Min"])*as.numeric(complex.rightoperand[i])
								complex.max.temp       <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Max"])*as.numeric(complex.rightoperand[i])
								result.i               <- matrix(c(complex.mat[i,"ParameterName"],complex.min.temp,complex.max.temp,"Unif*Constant"),nrow=1,ncol=4)
								colnames(result.i)     <- c("Parameter","Min","Max","DistClass")
								complex.dist.list[[i]] <-  as.data.frame(result.i)
							}
							### uniform simple parameter plus a constant
							if(complex.operators[i] == "+"){
								complex.min.temp       <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Min"])+as.numeric(complex.rightoperand[i])
								complex.max.temp       <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Max"])+as.numeric(complex.rightoperand[i])
								result.i               <- matrix(c(complex.mat[i,"ParameterName"],complex.min.temp,complex.max.temp,"Unif+Constant"),nrow=1,ncol=4)
								colnames(result.i)     <- c("Parameter","Min","Max","DistClass")
								complex.dist.list[[i]] <-  as.data.frame(result.i)
							}
							### uniform simple parameter minus a constant
							if(complex.operators[i] == "-"){
								complex.min.temp       <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Min"])-as.numeric(complex.rightoperand[i])
								complex.max.temp       <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Max"])-as.numeric(complex.rightoperand[i])
								result.i               <- matrix(c(complex.mat[i,"ParameterName"],complex.min.temp,complex.max.temp,"Unif-Constant"),nrow=1,ncol=4)
								colnames(result.i)     <- c("Parameter","Min","Max","DistClass")
								complex.dist.list[[i]] <-  as.data.frame(result.i)
							}
						}
					}
				} else {
					if(complex.leftoperand[i] %in% simple.pars & complex.rightoperand[i] %in% simple.pars){
						if(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Distribution"] == "logunif" | parameters.mat[grep(complex.rightoperand[i],simple.pars,fixed=T),"Distribution"] == "logunif"){
							complex.dist.list <- NULL
						} else {
							if(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Distribution"] == "lunif" | parameters.mat[grep(complex.rightoperand[i],simple.pars,fixed=T),"Distribution"] == "unif"){
								### uniform simple parameter divided by a uniform simple parameter
								if(complex.operators[i] == "/"){
									complex.min1.temp       <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Min"])
									complex.max1.temp       <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Max"])
									complex.min2.temp       <- as.numeric(parameters.mat[grep(complex.rightoperand[i],simple.pars,fixed=T),"Min"])
									complex.max2.temp       <- as.numeric(parameters.mat[grep(complex.rightoperand[i],simple.pars,fixed=T),"Max"])
									result.i               <- rbind(c(complex.mat[i,"ParameterName"],complex.min1.temp,complex.max1.temp,"Unif/Unif"),c(complex.mat[i,"ParameterName"],complex.min2.temp,complex.max2.temp,"Unif/Unif"))
									colnames(result.i)     <- c("Parameter","Min","Max","DistClass")
									complex.dist.list[[i]] <- as.data.frame(result.i)
								}
								### uniform simple parameter multiplied by a uniform simple parameter
								if(complex.operators[i] == "*"){
									complex.min1.temp       <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Min"])
									complex.max1.temp       <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Max"])
									complex.min2.temp       <- as.numeric(parameters.mat[grep(complex.rightoperand[i],simple.pars,fixed=T),"Min"])
									complex.max2.temp       <- as.numeric(parameters.mat[grep(complex.rightoperand[i],simple.pars,fixed=T),"Max"])
									result.i               <- rbind(c(complex.mat[i,"ParameterName"],complex.min1.temp,complex.max1.temp,"Unif*Unif"),c(complex.mat[i,"ParameterName"],complex.min2.temp,complex.max2.temp,"Unif*Unif"))
									colnames(result.i)     <- c("Parameter","Min","Max","DistClass")
									complex.dist.list[[i]] <- as.data.frame(result.i)
								}
								### uniform simple parameter plus a uniform simple parameter
								if(complex.operators[i] == "+"){
									complex.min1.temp       <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Min"])
									complex.max1.temp       <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Max"])
									complex.min2.temp       <- as.numeric(parameters.mat[grep(complex.rightoperand[i],simple.pars,fixed=T),"Min"])
									complex.max2.temp       <- as.numeric(parameters.mat[grep(complex.rightoperand[i],simple.pars,fixed=T),"Max"])
									result.i               <- rbind(c(complex.mat[i,"ParameterName"],complex.min1.temp,complex.max1.temp,"Unif+Unif"),c(complex.mat[i,"ParameterName"],complex.min2.temp,complex.max2.temp,"Unif+Unif"))
									colnames(result.i)     <- c("Parameter","Min","Max","DistClass")
									complex.dist.list[[i]] <- as.data.frame(result.i)
								}
								### uniform simple parameter minus a uniform simple parameter
								if(complex.operators[i] == "-"){
									complex.min1.temp      <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Min"])
									complex.max1.temp      <- as.numeric(parameters.mat[grep(complex.leftoperand[i],simple.pars,fixed=T),"Max"])
									complex.min2.temp      <- as.numeric(parameters.mat[grep(complex.rightoperand[i],simple.pars,fixed=T),"Min"])
									complex.max2.temp      <- as.numeric(parameters.mat[grep(complex.rightoperand[i],simple.pars,fixed=T),"Max"])
									result.i               <- rbind(c(complex.mat[i,"ParameterName"],complex.min1.temp,complex.max1.temp,"Unif-Unif"),c(complex.mat[i,"ParameterName"],complex.min2.temp,complex.max2.temp,"Unif-Unif"))
									colnames(result.i)     <- c("Parameter","Min","Max","DistClass")
									complex.dist.list[[i]] <- as.data.frame(result.i)
								}
							} else {
								### If none of the previous conditions were true set the ith entry to NULL
								complex.dist.list[[i]] <- NULL
							}
						}
					}
				}
			}
			if(any(lengths(complex.dist.list)>0)){
				complex.dist.df <- do.call(rbind,complex.dist.list)
				mode(complex.dist.df[,"Min"]) <- "numeric"
				mode(complex.dist.df[,"Max"]) <- "numeric"
				rownames(complex.dist.df)     <- NULL
			} else {
				complex.dist.df <- NULL
			}
		} else {
			complex.dist.list <- NULL
			complex.dist.df   <- NULL
		}
	} else {
		complex.dist.list <- NULL
		complex.dist.df   <- NULL
	}
	#print(complex.dist.df)
	##### Initial Trial Plot used for calculating space necessary for popnames and popsize labels.
	ymax.temp <- max(shapes.mat[,c("y_top_left","y_top_right")])*plotmargin[3]
	xmax.temp <- max(shapes.mat[,c("x_right")])*plotmargin[4]
	ymin.temp <- 0 + ymax.temp*(1-plotmargin[1])
	xmin.temp <- 0 + xmax.temp*(1-plotmargin[2])
	# Defining the blank plotting area
	par(mar=c(3,0.5,3,1.5),oma=c(2,2,2,2))
	plot(x=c(xmin.temp, xmax.temp), y=c(ymin.temp,ymax.temp), col="white",xlab="", ylab="",axes=F)
	### Adjusts plotmargin[2] to be sure that popnames and popsizes fit in plotting area.
	if(is.null(popnames)){
		### I changed this to require popnames, such that popnames = 1 and popnames = NULL do the same thing and will use "pop0", "pop1", etc.
		popnames.labels <- paste0("pop",c(0:(length(popsizes.present)-1)))
		xleft.adjust    <- max(strwidth(names(popsizes.present)))
	} else {
		if(popnames[1]==1 & length(popnames)==1){
			popnames.labels <- paste0("pop",c(0:(length(popsizes.present)-1)))
			xleft.adjust <- max(strwidth(names(popsizes.present)) + strwidth(popnames.labels))
		} else { 
			if(length(popnames)!=length(popsizes.present)){
				warning("number of popnames should equal the number of population samples on line 2 of template file")
				xleft.adjust <- max(strwidth(names(popsizes.present)))
			} else {
				popnames.labels <- popnames
				xleft.adjust <- max(strwidth(names(popsizes.present)) + strwidth(popnames.labels))
			}
		}
	}
	### limits for the plotting area
	ymax <- max(shapes.mat[,c("y_top_left","y_top_right")])*plotmargin[3]
	xmax <- max(shapes.mat[,c("x_right")])*plotmargin[4]
	ymin <- 0 + ymax*(1-plotmargin[1])
	xmin <- (0 - xleft.adjust) + xmax*(1-plotmargin[2])
	#### Do not move the next section before the test plotting section.
	if(show.priors & priors.panel.width > 0){
		### Setting the layout of the plotting area
		# Determines fraction of horizontal plotting area to use for plotting density functions of priors
		# number of graphical columns for the priors panel
		ncol.var     <- round(100*priors.panel.width)
		# number of graphical columns for the model
		ncol.model   <- (100-ncol.var)
		simple.par.per.type.list        <- lapply(X=list(popsizes.pars,time.pars,migrationrate.pars,migrants.pars,newsize.pars,growthrate.pars),FUN=function(x){intersect(x,simple.pars)})
		
		if(!is.null(complex.dist.df)){
			complex.par.per.type.list        <- lapply(X=list(popsizes.pars,time.pars,migrationrate.pars,migrants.pars,newsize.pars,growthrate.pars),FUN=function(x){intersect(x,complex.pars)})
			names(complex.par.per.type.list) <- c("popsizes","time","migrationrate","migrants.pars","newsize.pars","growthrate.pars")
			complex.pars.to.graph.temp       <- unique(complex.dist.df[,"Parameter"])
			complex.pars.to.graph            <- lapply(X=complex.par.per.type.list,FUN=function(x){intersect(x,complex.pars.to.graph.temp)})
			all.pars.to.graph                <- lapply(X=1:length(simple.par.per.type.list),FUN=function(x){union(complex.pars.to.graph[[x]],simple.par.per.type.list[[x]])})
		} else {
			all.pars.to.graph     <- simple.par.per.type.list
			complex.pars.to.graph <- NULL
		}
		other.simple.and.complex.list <- list(c(other.simple.pars,complex.pars))
		num.pars.per.type <- sapply(X=c(all.pars.to.graph,other.simple.and.complex.list),FUN=length)
		num.partypes      <- length(num.pars.per.type[num.pars.per.type!=0])
		numplots <- num.partypes
		# construct matrix defining layout of density plots
		layout.variables <- matrix(rep(2:(numplots+1),ncol.var),ncol=ncol.var)
		# construct matrix for the model part of the layout
		layout.model     <- matrix(data=1,ncol=ncol.model,nrow=nrow(layout.variables))
		# Define layout matrix
		layout.mat       <- cbind(layout.model,layout.variables)
		# Set plotting layout
		layout(layout.mat)
	}
	### For the case of all populations collapsed as a single population at all times, plot a single rectangle.
	if(all(shapes.mat[,"x_right"] %in% c(0,-1))){
		shapes.mat.temp <- unique(shapes.mat[,-1])
		plot(x=c(xmin, 10), y=c(ymin, ymax), col="white",xlab="", ylab="",axes=F)
		polygon(x=c(0,10,10,0), y=c(min(as.numeric(shapes.mat.temp[,"y_bottom_left"])),min(as.numeric(shapes.mat.temp[,"y_bottom_right"])),max(as.numeric(shapes.mat.temp[,"y_top_right"])),max(as.numeric(shapes.mat.temp[,"y_top_left"]))),col="gray",border="darkgray")
	} else {
		###
		# Defining the blank plotting area for model
		plot(x=c(xmin, xmax), y=c(ymin, ymax), col="white",xlab="", ylab="",axes=F)
		for(i in 1:nrow(shapes.mat)){
			shape.i <- shapes.mat[i,]
			if(shape.i["x_right"] %in% c(-1:0)){
				next
			}
			polygon(x=c(shape.i["x_left"],shape.i["x_right"],shape.i["x_right"],shape.i["x_left"]), y=c(shape.i["y_bottom_left"],shape.i["y_bottom_right"],shape.i["y_top_right"],shape.i["y_top_left"]),col="gray",border="darkgray")
		}
	}
	# ggplot2::ggplot() + ggplot2::scale_x_continuous(name="x") + ggplot2::scale_y_continuous(name="y") + ggplot2::geom_rect(data=test.df,mapping=ggplot2::aes(xmin= x_left,xmax= x_right,ymin= y_bottom_left,ymax= y_top_left))
	### Add arrows to show migration
	if(num.mig.matrices > 0){
		for(i in 1:length(period.sizes.v3)){
			period.mat     <- period.sizes.v3[[i]]
			shapes.mat.i   <- shapes.mat[c(shapes.mat[,"period"] %in% i),,drop=F]
			# skip to next i if only one population during the period
			if(nrow(shapes.mat.i)==1){
				next
			}
			if(period.mat["end.time","time"] %in% c(-1)){
				next
			}
			# x (time) coordinate where arrows will be drawn if there is migration during the period
			xmids.i   <- mean(range(shapes.mat.i[,c("x_left","x_right")]))
			xlefts.i  <- shapes.mat.i[,c("x_left")]
			xwidths.i <- diff(range(shapes.mat.i[,c("x_left","x_right")]))
			# ymids of each population (shape)
			ymids.i        <- sapply(X=1:nrow(shapes.mat.i),FUN=function(x){mean(range(shapes.mat.i[x,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right")]))})
			names(ymids.i) <- rownames(shapes.mat.i)
			### Get the migration matrix that applies to the region
			if(period.mat[1,"time"] %in% events.mat.values[,"time"]){
				# migration matrix index
				mig.mat          <- events.mat.values[which(events.mat.values[,"time"] == period.mat[1,"time"]),"migration matrix"]
				## There should never be more than one migration matrix in a time period.
				if(length(unique(mig.mat))>1){
					conincident.time.variables     <- events.mat[which(events.mat.values[,"time"] == period.mat[1,"time"]),"time"]
					conflicting.migration.matrices <- unique(mig.mat)
					stop(paste("Invalid model. Time variables:", paste(conincident.time.variables, collapse=" "),"have same value but different migration matrices:",paste(unique(mig.mat),collapse=" ")))
				}
				mig.mat          <- unique(mig.mat)
				# Corresponding migration matrix
				migration.matrix <- migration.matrices[[mig.mat+1]]
			} else {
				# The first migration matrix is the one used for time 0, if it exists, unless the events matrix specifies otherwise.
				migration.matrix <- migration.matrices[[1]]
			}
			# if migration matrix is not the zero matrix, get row and column indices of the nonzero values. Probably should use names instead of numbers.
			if(!all(migration.matrix==0)){
				ymids.source.mat.i      <- matrix(data=rep(ymids.i,length(ymids.i)), nrow=length(ymids.i),ncol=length(ymids.i),byrow=F)
				ymids.sink.mat.i        <- matrix(data=rep(ymids.i,length(ymids.i)), nrow=length(ymids.i),ncol=length(ymids.i),byrow=T)
				source.sink.mat         <- cbind(pop.source=rep( rownames(shapes.mat.i), each=length(ymids.i)), pop.sink=rep(rownames(shapes.mat.i),length(ymids.i)),ymids.source = c(ymids.source.mat.i), ymids.sink = c(ymids.sink.mat.i), migration.parameter = c(migration.matrix[rownames(shapes.mat.i),rownames(shapes.mat.i)]))
				source.sink.mat.ordered <- source.sink.mat[order(source.sink.mat[,"pop.source"], source.sink.mat[,"pop.sink"]),]
				pop.pair.mig            <- do.call(rbind,lapply(X=1:nrow(source.sink.mat.ordered),FUN=function(x){paste(paste(unname(sort(source.sink.mat[x,c("pop.source","pop.sink")])),collapse=" "),source.sink.mat[x,c("migration.parameter")],collapse=" ")}))
				source.sink.df          <- as.data.frame(source.sink.mat.ordered)
				mode(source.sink.df[,"ymids.source"]) <- "numeric"
				mode(source.sink.df[,"ymids.sink"])   <- "numeric"
				source.sink.df[,"poppair.mig"]        <- pop.pair.mig
				if(all(source.sink.df$migration.parameter==0)){
					next
				}
				source.sink.df                <- source.sink.df[which(source.sink.df$migration.parameter!=0),]
				x.adj.tol                     <- 0.15
				x.adj                         <- seq(from=(-x.adj.tol),to=x.adj.tol ,length.out=length(unique(source.sink.df[,"poppair.mig"])))
				source.sink.df[,"xadj"]       <- x.adj[match(source.sink.df[,"poppair.mig"],unique(source.sink.df[,"poppair.mig"]))]
				mode(source.sink.df[,"xadj"]) <- "numeric"
				source.sink.df[,"x"]          <- (xmids.i+(max(shapes.mat.i[,"x_right"]-shapes.mat.i[,"x_left"])*source.sink.df[,"xadj"]))
				apply(X=source.sink.df[,c("ymids.source","ymids.sink")],MARGIN=1,FUN=mean)
				source.sink.df[,"ytext"] <- apply(X=source.sink.df[,c("ymids.source","ymids.sink")],MARGIN=1,FUN=mean)
				arrows(x0=source.sink.df[,"x"],y0=source.sink.df[,"ymids.source"],y1=source.sink.df[,"ymids.sink"],length=0.05,col="black")
				text.df <- unique(source.sink.df[,c("migration.parameter","x","ytext")])
				text(x=text.df[,"x"],y=text.df[,"ytext"],labels=text.df[,"migration.parameter"],adj=c(0.5,-0.5),cex=label.cex,srt=90)
			}
		} ### end arrows loop
	} ### end if statement
	### Add dashed vertical lines to separate time periods at historical events
	segments(x0=events.mat.values[,"time"],y0=ymin,y1=ymax,lty=2)
	### Add dashed vertical line at t=0
	segments(x0=0,y0=ymin,y1=ymax,lty=2)
	### Add text labels for event times except for events that occur at time zero...
	if(num.events>0){
		if(any(events.mat.values[,"time"] != 0)){
			events.mark     <- which(events.mat.values[,"time"] != 0)
			xylab.event.mat <- cbind(x=events.mat.values[events.mark,"time"],y=rep(mean(c(0,ymin)),length(events.mark)),lab=events.mat[events.mark,"time"])
			### identify and merge labels for events that occur at the same time
			while(any(table(xylab.event.mat[,"x"])>1)){
				num.concurrent.events <- length(which(table(xylab.event.mat[,"x"])>1))
				time.events     <- names(which(table(xylab.event.mat[,"x"])>1)[1])
				merged.lab      <- paste0(xylab.event.mat[which(xylab.event.mat[,"x"]==time.events),"lab"],collapse=", ")
				merged.row      <- cbind(unique(xylab.event.mat[which(xylab.event.mat[,"x"]==time.events),c("x","y"),drop=F]),merged.lab)
				xylab.event.mat <- rbind(xylab.event.mat[-which(xylab.event.mat[,"x"]==time.events),],merged.row)
			}
			if(nrow(xylab.event.mat)>1){
				xylab.event.mat[which((1:nrow(xylab.event.mat) / 2 ) == round((1:nrow(xylab.event.mat) / 2 ))),"y"] <- mean(c(0,0,0,ymin))
			}
			xylab.event.mat[which((1:nrow(xylab.event.mat) / 2 ) != round((1:nrow(xylab.event.mat) / 2 ))),"y"] <- mean(c(0,ymin,ymin,ymin))
			### plot text labels for event times. Even and odd rows are plotted separately
			text(x=as.numeric(xylab.event.mat[,"x"]), y=as.numeric(xylab.event.mat[,"y"]), labels=paste0(" ",xylab.event.mat[,"lab"]), adj= c(0,0.5), cex=label.cex)
		} else {
			xylab.event.mat <- NULL
		}
	} else {
		xylab.event.mat <- NULL
	}
	### Add time 0 label. Labels of time parameters that equal zero should be merged with the "0" label
	if(any(events.mat.values[,"time"] == 0)){
		events.mark.0     <- which(events.mat.values[,"time"] == 0)
		xylab.event.0.mat <- cbind(x=events.mat.values[events.mark.0,"time"],y=rep(mean(c(0,ymin)),length(events.mark.0)),lab=events.mat[events.mark.0,"time"])
		xylab.event.0.mat <- unique(rbind(c(x=0,y=mean(c(0,ymin)),lab="0"),xylab.event.0.mat))
		while(any(table(xylab.event.0.mat[,"x"])>1)){
			merged.lab      <- paste0(xylab.event.0.mat[,"lab"],collapse=", ")
			xylab.event.0.mat <- cbind(x=0,y=mean(c(0,ymin)),lab=merged.lab)
		}
		### If the merged time zero label is still just "0", replace it with " 0"
		if(xylab.event.0.mat[,"lab"]=="0"){
			xylab.event.0.mat[,"lab"] <- " 0"
		}
		### plot text label for time 0 and for events that occur at time 0
		text(x=as.numeric(xylab.event.0.mat[,"x"]), y=as.numeric(xylab.event.0.mat[,"y"]), labels=xylab.event.0.mat[,"lab"], adj= c(-0.1,0.5), cex=label.cex)
	} else {
		text(x=0,y=mean(c(0,ymin)),labels=" 0",adj= c(-0.1,0.5),cex=label.cex)
		xylab.event.0.mat <- NULL
	}
	### Add labels for initial population sizes
	shapes.mat1 <- shapes.mat[shapes.mat[,"period"] == 1,]
	y.poplabels <- sapply(X=1:nrow(shapes.mat1),FUN=function(x){mean(range(shapes.mat1[x,c("y_top_left","y_top_right","y_bottom_left","y_bottom_right")]))})
	text(x=0,y=y.poplabels,labels=names(popsizes.present),adj= c(1.2,0.5),cex=label.cex)
	### Add labels for popnames
	if(!is.null(popnames)){
		text(x=xmin,y=y.poplabels,labels=popnames.labels,adj= c(0,0.5),cex=label.cex)
	}
	### Add labels for growthrate change parameters, migration events, and newsize events
	if(length(growthrate.pars)>0){
		times.at.growthrate <- events.mat.values[grep(growthrate.pars,events.mat[,"growthrate"],fixed=T),"time"]
	} else {
		times.at.growthrate <- NULL
	}
	if(length(migrants.pars)>0){
		times.at.migrants <- events.mat.values[grep(migrants.pars,events.mat[,"migrants"],fixed=T),"time"]
	} else {
		times.at.migrants <- NULL
	}
	if(length(newsize.pars)>0){
		times.at.newsize <- events.mat.values[grep(newsize.pars,events.mat[,"newsize"],fixed=T),"time"]
	} else {
		times.at.newsize <- NULL
	}
	time.at.event.mat    <- cbind(time=c(times.at.growthrate,times.at.migrants,times.at.newsize),parameter=c(growthrate.pars,migrants.pars,newsize.pars))
	xylab.event.mat.all  <- rbind(xylab.event.0.mat, xylab.event.mat)
	if(!is.null(time.at.event.mat) & !is.null(xylab.event.mat.all)){
		rownames(xylab.event.mat.all) <- NULL
		for(i in 1:nrow(time.at.event.mat)){
			event.i            <- grep(time.at.event.mat[i,"time"],xylab.event.mat.all[,"x"])
			timelabel.i.height <- strheight(xylab.event.mat.all[event.i,"lab"])
			event.i.y          <- as.numeric(xylab.event.mat.all[event.i,"y"])-timelabel.i.height
			event.i.x          <- as.numeric(xylab.event.mat.all[event.i,"x"])
			label.i            <- time.at.event.mat[i,"parameter"]
			text(x=event.i.x, y=event.i.y, labels=paste0(" ",label.i), adj= c(0,1), cex=label.cex)
			xylab.event.mat.all <- rbind(xylab.event.mat.all,c(event.i.x,event.i.y,label.i))
		}
	}
	### Add x axis label
	mtext(text="generations ago (increasing -->) ",side=1)
	### Add title
	box(which="figure")
	mtext(text=model.title,outer=T,line=0.3,adj=c(mean(c(0,1-(priors.panel.width)))))
	
	###### Experimenting here #####
	if(show.priors & priors.panel.width > 0){
		## Determine max and min values with nonzero probability for complex parameters that should be plotted
		if(!is.null(complex.dist.df)){
			if(any(complex.dist.df[,"DistClass"] %in% c("Unif/Constant","Unif*Constant","Unif+Constant","Unif-Constant"))){
				#complex.ranges.df.temp <- complex.dist.df[complex.dist.df[,"DistClass"] %in% c("Unif/Constant","Unif*Constant","Unif+Constant","Unif-Constant"),]
				complex.ranges.df           <- complex.dist.df[complex.dist.df[,"DistClass"] %in% c("Unif/Constant","Unif*Constant","Unif+Constant","Unif-Constant"),]
				colnames(complex.ranges.df) <- c("Parameter","Min","Max","DistClass")
				mode(complex.ranges.df[,"Min"]) <- "numeric"
				mode(complex.ranges.df[,"Max"]) <- "numeric"
			} else {
				complex.ranges.df <- NULL
			}
			if(any(complex.dist.df[,"DistClass"] %in% "Unif*Unif")){
				complex.ranges.df.temp   <- complex.dist.df[complex.dist.df[,"DistClass"] %in% "Unif*Unif",]
				colnames(complex.ranges.df.temp) <- c("Parameter","Min","Max","DistClass")
				mode(complex.ranges.df.temp[,"Min"]) <- "numeric"
				mode(complex.ranges.df.temp[,"Max"]) <- "numeric"
				complex.pars.temp        <- unique(complex.ranges.df.temp[,"Parameter"])
				for(i in 1:length(complex.pars.temp)){
					param.temp <- complex.pars.temp[i]
					rows.temp  <- which(complex.ranges.df.temp[,"Parameter"]==param.temp)
					min.temp   <- prod(complex.ranges.df.temp[rows.temp,"Min"])
					max.temp   <- prod(complex.ranges.df.temp[rows.temp,"Max"])
					#complex.ranges.df <- rbind(complex.ranges.df, c(param.temp,min.temp,max.temp,"Unif*Unif"))
					complex.ranges.df <- rbind(complex.ranges.df, data.frame('Parameter'=param.temp,'Min'=min.temp,'Max'=max.temp,'DistClass'="Unif*Unif"))
				}
				#colnames(complex.ranges.df) <- c("Parameter","Min","Max","DistClass")
				#mode(complex.ranges.df[,"Min"]) <- "numeric"
				#mode(complex.ranges.df[,"Max"]) <- "numeric"
			}
			if(any(complex.dist.df[,"DistClass"] %in% "Unif+Unif")){
				complex.ranges.df.temp <- complex.dist.df[complex.dist.df[,"DistClass"] %in% "Unif+Unif",]
				colnames(complex.ranges.df.temp) <- c("Parameter","Min","Max","DistClass")
				mode(complex.ranges.df.temp[,"Min"]) <- "numeric"
				mode(complex.ranges.df.temp[,"Max"]) <- "numeric"
				complex.pars.temp       <- unique(complex.ranges.df.temp[,"Parameter"])
				for(i in 1:length(complex.pars.temp)){
					param.temp        <- complex.pars.temp[i]
					rows.temp         <- which(complex.ranges.df.temp[,"Parameter"]==param.temp)
					min.temp          <- sum(complex.ranges.df.temp[rows.temp,"Min"])
					max.temp          <- sum(complex.ranges.df.temp[rows.temp,"Max"])
					#complex.ranges.df <- rbind(complex.ranges.df,c(param.temp,min.temp,max.temp,"Unif+Unif"))
					complex.ranges.df <- rbind(complex.ranges.df, data.frame('Parameter'=param.temp,'Min'=min.temp,'Max'=max.temp,'DistClass'="Unif+Unif"))
				}
				#colnames(complex.ranges.df) <- c("Parameter","Min","Max","DistClass")
				#mode(complex.ranges.df[,"Min"]) <- "numeric"
				#mode(complex.ranges.df[,"Max"]) <- "numeric"
			}
		} else {
			complex.ranges.df <- NULL
		}
		### Data frame with range, distribution type, and parameter type for both simple and complex parameters to be graphed
		if(!is.null(parameters.mat)){
			simple.ranges.df <- as.data.frame(parameters.mat[,c("ParameterName","Min","Max","Distribution")])
			colnames(simple.ranges.df) <- c("Parameter","Min","Max","DistClass")
			mode(simple.ranges.df[,"Min"]) <- "numeric"
			mode(simple.ranges.df[,"Max"]) <- "numeric"
			if(!is.null(complex.ranges.df)){
				ranges.df <- rbind(simple.ranges.df,complex.ranges.df)
			} else {
				ranges.df <- simple.ranges.df
			}
			ranges.df[,"ParameterType"] <- ""
			if(any(ranges.df[,"Parameter"] %in% time.pars)){
				ranges.df[ranges.df[,"Parameter"] %in% time.pars,"ParameterType"] <- "time"
				min.time <- min(ranges.df$Min[ranges.df$ParameterType %in% "time"])
				max.time <- min(ranges.df$Max[ranges.df$ParameterType %in% "time"])
			}
			if(any(ranges.df[,"Parameter"] %in% popsizes.pars)){
				ranges.df[ranges.df[,"Parameter"] %in% popsizes.pars,"ParameterType"] <- "popsizes"
				min.popsizes <- min(ranges.df$Min[ranges.df$ParameterType %in% "popsizes"])
				max.popsizes <- min(ranges.df$Max[ranges.df$ParameterType %in% "popsizes"])
			}
			if(any(ranges.df[,"Parameter"] %in% migrationrate.pars)){
				ranges.df[ranges.df[,"Parameter"] %in% migrationrate.pars,"ParameterType"] <- "migrationrate"
					min.migrationrate <- min(ranges.df$Min[ranges.df$ParameterType %in% "migrationrate"])
				max.migrationrate <- min(ranges.df$Max[ranges.df$ParameterType %in% "migrationrate"])
			}
			if(any(ranges.df[,"Parameter"] %in% migrants.pars)){
				ranges.df[ranges.df[,"Parameter"] %in% migrants.pars,"ParameterType"] <- "migrants"
				min.migrants <- min(ranges.df$Min[ranges.df$ParameterType %in% "migrants"])
				max.migrants <- min(ranges.df$Max[ranges.df$ParameterType %in% "migrants"])
			}
			if(any(ranges.df[,"Parameter"] %in% newsize.pars)){
				ranges.df[ranges.df[,"Parameter"] %in% newsize.pars,"ParameterType"] <- "newsize"
				min.newsize <- min(ranges.df$Min[ranges.df$ParameterType %in% "newsize"])
				max.newsize <- min(ranges.df$Max[ranges.df$ParameterType %in% "newsize"])
	
			}
			if(any(ranges.df[,"Parameter"] %in% growthrate.pars)){
				ranges.df[ranges.df[,"Parameter"] %in% growthrate.pars,"ParameterType"] <- "growthrate"
				min.growthrate <- min(ranges.df$Min[ranges.df$ParameterType %in% "growthrate"])
				max.growthrate <- min(ranges.df$Max[ranges.df$ParameterType %in% "growthrate"])

			}
			if(any(ranges.df[,"Parameter"] %in% other.pars)){
				ranges.df[ranges.df[,"Parameter"] %in% other.pars,"ParameterType"] <- "other"
				min.other <- min(ranges.df$Min[ranges.df$ParameterType %in% "other"])
				max.other <- min(ranges.df$Max[ranges.df$ParameterType %in% "other"])
			}
		} else {
			ranges.df <- NULL
		}
		#######################################
		####### Probability Distributions #####
		#######################################
		###### Better way to get pdf vals
		if(!is.null(ranges.df)){
			dist.df <- as.data.frame(matrix(data=0,ncol=nrow(ranges.df),nrow=500))
			colnames(dist.df) <- ranges.df[,"Parameter"]
			rownames(dist.df) <- NULL
			#mode(dist.df) <- "numeric"
			for(i in 1:nrow(ranges.df)){
				type.temp <- ranges.df[i,"ParameterType"]
				## Min and Max among parameters of the same type
				min.temp  <- min(ranges.df[ranges.df[,"ParameterType"] %in% type.temp,"Min"])
				max.temp  <- max(ranges.df[ranges.df[,"ParameterType"] %in% type.temp,"Max"])
				x0        <- min.temp-abs((min.temp*0.1))
				x1        <- max.temp+abs((max.temp*0.1))
				z.samples <- seq(from=x0,to=x1,length.out=500)
				if(ranges.df[i,"DistClass"] %in% c("unif","Unif*Constant","Unif/Constant","Unif+Constant","Unif-Constant")){
					#dist.df.col.i <- dunif(z.samples,min=ranges.df[i,"Min"],max=ranges.df[i,"Max"])
					#dist.df[,i]   <- dist.df.col.i
					dist.df[,i]   <- dunif(z.samples,min=ranges.df[i,"Min"],max=ranges.df[i,"Max"])
				}
				if(ranges.df[i,"DistClass"] %in% c("logunif","Logunif*Constant","Logunif/Constant","Logunif+Constant","Logunif-Constant")){
					#dist.df.col.i <- dunif(z.samples,min=ranges.df[i,"Min"],max=ranges.df[i,"Max"])
					#dist.df[,i]   <- dist.df.col.i
					dist.df[,i]   <- dlogunif(z.samples,min=ranges.df[i,"Min"],max=ranges.df[i,"Max"])
				}
				if(ranges.df[i,"DistClass"] == "Unif+Unif"){
					a.vals      <- complex.dist.df[complex.dist.df[,"Parameter"] %in% ranges.df[i,"Parameter"],"Min"]
					b.vals      <- complex.dist.df[complex.dist.df[,"Parameter"] %in% ranges.df[i,"Parameter"],"Max"]
					dist.df[,i] <- sapply(z.samples,FUN=function(x){sumunif(z=x,a=a.vals,b=b.vals)})
				}
				if(ranges.df[i,"DistClass"] == "Unif*Unif"){
					a.vals      <- complex.dist.df[complex.dist.df[,"Parameter"] %in% ranges.df[i,"Parameter"],"Min"]
					b.vals      <- complex.dist.df[complex.dist.df[,"Parameter"] %in% ranges.df[i,"Parameter"],"Max"]
					dist.df[,i] <- sapply(z.samples,FUN=function(x){produnif(z=x,a=a.vals,b=b.vals)})
				}
				if(ranges.df[i,"DistClass"] %in% c("Unif-Unif","Unif/Unif","Logunif*Unif","Logunif/Unif","Logunif+Unif","Logunif-Unif","Logunif*Logunif","Logunif/Logunif","Logunif+Logunif","Logunif-Logunif")){
					next
				}
			}
		}
		############
		### Graphing the Probability Density Curves
		############
		linecols    <- c("red","blue","darkgreen","darkorange","purple","black","gray","yellow","red")
		line.types  <- c(1,1,3,4,1,1,1,1,3,1,1,1)
		line.widths <- c(2,1,1,1,0.75,0.75,0.75,0.75,0.75,0.75)
		## Extra space to add above density plots for including variable lables. Value of 1 leaves no space.
		varnames.spacer <- 2
		panel.title <- "Parameter Probability Densities"
		#pars.to.plot.list <- list(popsizes.pars.temp,time.pars.temp,migrationrate.pars.temp,migrants.pars.temp,newsize.pars.temp,growthrate.pars.temp)
		pars.to.plot.list <- c("popsizes","time","migrationrate","migrants","newsize","growthrate")
		xaxis.label <- c("effective pop. sizes (t = 0)","generations ago","migration rate","migrants (fraction of source population)","population size change (new/old)","population growth rate")
		for(z in 1:length(pars.to.plot.list)){
			if(z==1){
				par(mar=c(3,2,0,1))
			} else {
				par(mar=c(3,2,0,1))
			}
		#	pars.temp <- pars.to.plot.list[[z]]
			if(any(ranges.df[,"ParameterType"] %in% pars.to.plot.list[[z]])){
				pars.temp <- ranges.df[(ranges.df[,"ParameterType"] %in% pars.to.plot.list[[z]]),"Parameter"]
				pars.min  <- min(ranges.df[(ranges.df[,"ParameterType"] %in% pars.to.plot.list[[z]]),"Min"])
				pars.max  <- max(ranges.df[(ranges.df[,"ParameterType"] %in% pars.to.plot.list[[z]]),"Max"])
				x0        <- pars.min-abs((pars.min*0.1))
				x1        <- pars.max+abs((pars.max*0.1))
				dist.df.temp   <- dist.df[,pars.temp,drop=F]
				pd.Max         <- max(dist.df.temp)
				if(pd.Max==0){
					next
				}
				ycurves_top    <- pd.Max
				ylabelarea_top <- ycurves_top*varnames.spacer
				z.samples <- seq(from=x0,to=x1,length.out=500)
				plot(x=c(x0,x1),y=c(0,ylabelarea_top),main="",xlab="",ylab="",axes=F,type="l",lty=1, col="white")
				box(which="figure")
				for(i in 1:ncol(dist.df.temp)){
					if(max(dist.df.temp[,i])==0){
						next
					}
					lines(z.samples, dist.df.temp[,i], col=linecols[i], lty=line.types[i],lwd=line.widths[i])
				}
				if(z==1){
					mtext(side=3,text=panel.title,font=1,cex=0.75,line=0.3)
				}
				xmin.label <- formatC(x0,format="g",digits=2)
				xmax.label <- formatC(x1,format="g",digits=2)
				axis(side=1,at=c(x0,x1),labels=c(xmin.label,xmax.label),cex.axis=0.8,hadj=(0.5),padj=(-1))
				mtext(text=xaxis.label[z],side=1,line=1.2,cex=0.6)
				plot_area_width  <- par("usr")[2]-par("usr")[1]
				plot_area_height <- par("usr")[4]-par("usr")[3]
				plot_area_ymax   <- par("usr")[4]
				plot_area_ymin   <- par("usr")[3]
				plot_area_xmax   <- par("usr")[2]
				plot_area_xmin   <- par("usr")[1]
				cex.temp <- 0.8
				label_heights           <- max(strheight(pars.temp,cex=cex.temp))
				label_widths            <- strwidth(pars.temp,cex=cex.temp)
				cumulative_label_widths <- cumsum(label_widths)
				### Place nth label to the left of the nth cumalative width
				ybottom_varlabels       <- ycurves_top*1.2
				x_varlabels             <- (cumulative_label_widths*1.2)+x0
				x_varlabels_right.edge  <- (max(cumulative_label_widths)*1.25)+x1
				# line color and type legend
				y_linelabels  <- (ybottom_varlabels + label_heights)*1.1
				x0_linelabels <- (x_varlabels-label_widths)
				x1_linelabels <- x_varlabels
				counter=0
				segments(x0=x0_linelabels,x1=x1_linelabels,y0=y_linelabels,col=linecols[1:length(pars.temp)],lty=line.types[1:length(pars.temp)],lwd=line.widths[i])
				text(x=x_varlabels,y=ybottom_varlabels,labels=pars.temp,col=linecols[1:length(pars.temp)],cex=cex.temp,adj=c(1,0))
			}
		}
		### Rather than plotting the other parameters, just show them in a text table
		other.and.complex.pars <- unique(c(other.pars,complex.pars))
		if(length(rules.mat) > 0){
			numlines <- length(other.and.complex.pars) + nrow(rules.mat) + 2
		} else {
			numlines <- length(other.and.complex.pars) + 1
		}
		if(length(other.and.complex.pars) > 0){
			par(mar=c(0,0,0,0))
			plot(x=0:numlines,y=0:numlines,col="white",axes=F,xlab="",ylab="",main="",cex.main=1)
			box(which="figure")
			# names of parameters in column 1
			if(any(parameters.mat[,"ParameterName"] %in% other.pars)){
				others.simple.mat <- parameters.mat[parameters.mat[,"ParameterName"] %in% other.pars,c("ParameterName","Distribution","Min","Max"),drop=F]
				others.simple.string <- list(); length(others.simple.string) <- nrow(others.simple.mat)
				for(i in 1:nrow(others.simple.mat)){
					others.simple.string[[i]] <- paste0(others.simple.mat[i,1]," ~ ",others.simple.mat[i,2],"(",others.simple.mat[i,3],",",others.simple.mat[i,4],")")
				}
				others.simple.string <- unlist(others.simple.string)
			} else {
				others.simple.string <- NULL
			}
			if(length(complex.pars)>0){
				complex.text.string.mat <- cbind(complex.mat[,c("ParameterName","function"),drop=F],complex.text.string)
				complex.string          <- unlist(apply(X=complex.text.string.mat,MARGIN=1,FUN=paste,collapse=" "))
			} else {
				complex.string <- NULL
			}
			if(length(rules.mat)>0){
				rules.string <- c("Rules",unlist(apply(X=rules.mat,MARGIN=1,FUN=paste,collapse=" ")))
				rules_height <- strheight(rules.string)
				rules_width  <- strwidth(rules.string)
			} else {
				rules.string <- NULL
			}
			others.complex.string <- c("Complex and other parameters",others.simple.string,complex.string)
			others.complex_height <- strheight(others.complex.string)
			others.complex_widths <- strwidth(others.complex.string)
			cex.temp <- 1
			y_others.complex <- (par("usr")[4]-(cumsum(others.complex_height)*1.4)) + others.complex_height[1]
			if(length(others.complex.string)>0){
				text(x=rep(0,length(others.complex.string)),y=(y_others.complex-0.5),labels=others.complex.string,cex=cex.temp,pos=4,font=c(2,rep(1,length(others.complex.string)-1)))
			}
			if(length(rules.string)>0){
				y_rules <- (min(y_others.complex)-(cumsum(rules_height)*1.4))-rules_height[1]
				text(x=0,y=(y_rules-0.5),labels=rules.string,cex=cex.temp,pos=4,font=c(2,rep(1,length(rules.string)-1)))
			}
		}
	}
	### hold the plot as an object
	result <- recordPlot()
	### return graphical settings back to their original state
	layout(1)
	### output the plot
	result
}
#' @examples
#' source("plot_model.R")
#' dev.new(width=10,height=6)
#' Model5.plot  <- plot_model(tpl.path="example_5.tpl",  est.path="example_5.est")
#' Model6.plot  <- plot_model(tpl.path="example_6.tpl",  est.path="example_6.est")
#' Model7.plot  <- plot_model(tpl.path="example_7.tpl",  est.path="example_7.est")
#' Model8.plot  <- plot_model(tpl.path="example_8.tpl",  est.path="example_8.est")
#' Model9.plot  <- plot_model(tpl.path="example_9.tpl",  est.path="example_9.est")
#' Model10.plot <- plot_model(tpl.path="example_10.tpl", est.path="example_10.est")
#' Model11.plot <- plot_model(tpl.path="example_11.tpl", est.path="example_11.est")
#' Model12.plot <- plot_model(tpl.path="example_12.tpl", est.path="example_12.est")
#' Model13.plot <- plot_model(tpl.path="example_13.tpl", est.path="example_13.est")
#' Model14.plot <- plot_model(tpl.path="example_14.tpl", est.path="example_14.est")
#' Model15.plot <- plot_model(tpl.path="example_15.tpl", est.path="example_15.est")
#' Model16.plot <- plot_model(tpl.path="example_16.tpl", est.path="example_16.est")
#' Model17.plot <- plot_model(tpl.path="example_17.tpl", est.path="example_17.est")
#' Model18.plot <- plot_model(tpl.path="example_18.tpl", est.path="example_18.est")
#' 
#' list.of.models <- list(Model5.plot,Model6.plot,Model7.plot,Model8.plot,Model9.plot,Model10.plot,Model11.plot,Model12.plot,Model13.plot,Model14.plot,Model15.plot,Model16.plot,Model17.plot,Model18.plot)
#' pdf("cartoon_models.pdf",width=10,height=6)
#' list.of.models
#' dev.off()



#' @title dlunif function
#' 
#' loguniform probability density function
#' 
#' @param x sample quantile
#' @param min minimum of sampling distribution
#' @param max maximum of sampling distribution
#' @param base base of exponent. Default exp(1)
#' @return probability density at x
#' @export dlogunif
dlogunif <- function (x, min, max, base = exp(1)) {
	return((1/(log(max, base) - log(min, base))) * (1/x))
}

#' @title Sample from loguniform distribution
#' 
#' 
#' @param n Number of samples to draw
#' @param min minimum of sampling distribution
#' @param max maximum of sampling distribution
#' @param base base of exponential parameter. Default exp(1).
rlogunif <- function (n, min, max, base = exp(1)){
	return(base^(runif(n, log(min, base),log(max, base))))
}

#' @title Probability density function for the quotient of two uniform distributions that both have minimum value at zero.
#' 
#' Probability density for z in Z ~ X1/X2, where X1 ~ U(0,b1) and X2 ~ U(0,b2), b1 > 0, b2 > 0.
#' 
#' @param z quantile
#' @param b Numerical vector with b1 and b2 for X1~Unif(0,b1) and X2~Unif(0,b2).
#' @return Probability density of z in Z
#' @export dunif.quotient.a0
dunif.quotient.a0 <- function(z,b){
	d.standard <- (0.5*((((sign(z-1)+1)/2)/(z^2))+((sign(1-z)+1)/2)))
	result     <- d.standard*(b[1]/b[2])
	result
}

#' @title Probability density function for the quotient of two identically independent uniform distributions.
#' 
#' Probability density for z in Z ~ X1/X2, where X1 & X2 ~ U(a,b), a & b in Reals, and a < b.
#' The function dunif.quotient.a0 does not require X1 and X2 to be identical, but has the requirement that the lower limit of both distributions equal 0.
#' 
#' @param z quantile
#' @param a Number indicating the lower limit of X1 and X2
#' @param b Number indicating the lower upper of X1 and X2
#' @return Probability density of z in Z
#' @export dunif.quotient.iid
dunif.quotient.iid <- function(z,a,b){
	if(a > b){
		stop("'a' must be < 'b'")
	}
	if(0 < a & a < b){
		if(z > (b/a) | z < (a/b)){
			result <- 0
		} else {
			if(z <= 1){
				result <- (((z^2)*(b^2))-(a^2)) / (2*(z^2)*((b-a)^2))
			}
			if(z > 1){
				result <- (((b^2)-((a^2)*(z^2)))/(2 * (z^2) * ((b-a)^2)))
			}
		}
	}
	if(a < 0 & 0 < b & abs(b) > abs(a)){
		if(z <= (b/a)){
			result <- ((b^2) + (a^2)) / (2 * (z^2) * ((b-a)^2))
		}
		if((b/a) < z & z <= (a/b)){
			result <- ((a^2)*((z^2) + 1)) / (2 * (z^2) * ((b-a)^2))
		}
		if((a/b) < z & z <= 1){
			result <- ((a^2) + (b^2)) / (2 * ((b-a)^2))
		}
		if(1 < z){
			result <- ((a^2) + (b^2)) / (2 * (z^2) * ((b-a)^2))
		}
	}
	if(a < 0 & 0 < b & abs(b) < abs(a)){
		if(z <= (a/b)){
			result <- ((a^2) + (b^2)) / (2 * (z^2) * ((b-a)^2))
		}
		if((a/b) < z & z <= (b/a)){
			((b^2) + (1 + (z^2))) / (2 * (z^2) * ((b-a)^2))
		}
		if((b/a) < z & z <= 1){
			result <- ((a^2) + (b^2)) / (2 * ((b-a)^2))
		}
		if(1 < z){
			result <- ((a^2) + (b^2)) / (2 * (z^2) * ((b-a)^2))
		}
	}
	if(a < b & b < 0){
		if(z > (a/b) | z < (b/a)){
			result <- 0
		} else {
			if(1 < z){
				result <- abs((((z^2)*(b^2))-(a^2)) / (2*(z^2)*((b-a)^2)))
			}
			if(z <= 1){
				result <- abs(((b^2)-((a^2)*(z^2))) / (2 * (z^2) * ((b-a)^2)))
			}
		}
	}
	result
}


#' @title PDF or CDF for the sum of n independent uniformly continuous random variables.
#' 
#' Evaluates the probability density function or cumulative distribution of Z=z for Z(z) = Xi+...+Xn, for i = 1...n, Xi ~ U(ai,bi), where 0 < ai < bi
#' The PDF and CDF is calculated using Equations 2.4 and 2.7, respectively, from Ishihara, Tatsuo. 2002. The Distribution of the Sum and the Product of Independent Uniform Random Variables Distributed at Different Intervals. Transactions - Japan Society for Industrial and Applied Mathematics, 12(3), 197-208.
#' 
#' @param z Number or numerical vector of quantiles
#' @param a Vector holding the mins of each X
#' @param b Vector holding the max of each X
#' @param cumulative Whether or not return the cumulative probability instead of the probability. Default is FALSE.
#' @return Number or numerical vector with either the probability density or the cumulative density for each input value z.
#' @export sumunif
sumunif <- function(z,a,b,cumulative=F){
	if(!cumulative & (z <= sum(a) | z > sum(b))){
		return(0)
	}
	if(cumulative & z <= sum(a)){
		return(0)
	}
	if(cumulative & z >= sum(b)){
		return(1)
	}
	n = length(a)
	powerSet.list               <- pset(1:n)
	solution.mat                <- matrix(data=NA,nrow=length(powerSet.list),ncol=9)
	colnames(solution.mat)      <- c("No.","alpha","(-1)Nalpha","Sn,alpha","z","(z-Sn,alpha)","((z-Sn,alpha)+)","((z-Sn,alpha)^(n-1))","((-1)^Nalpha)*((z-Sn,alpha)^(n-1))")
	solution.mat[,"No."]        <- 1:nrow(solution.mat)
	solution.mat[,"alpha"]      <- sapply(X = powerSet.list, FUN=paste, collapse = " ")
	solution.mat[,"(-1)Nalpha"] <- (-1)^lengths(powerSet.list)
	solution.mat[,"z"]          <- rep(z,nrow(solution.mat))
	C  <- b-a
	CN <- max(cumprod(C))
	if(cumulative){
		part1 <- 1/((factorial(n))*CN)
	} else {
		part1 <- 1/((factorial(n-1))*CN)
	}
	for(i in 1:nrow(solution.mat)){
		alpha.j <- rep(0,n)
		if(i>1){
			alpha.j[1:n %in% powerSet.list[[i]]] <- 1
		}
		Sn.i    <- sum(a + (alpha.j*C))
		solution.mat[i,"Sn,alpha"] <- Sn.i
	}
	solution.df        <- as.data.frame(solution.mat)
	mode(solution.df[,"(-1)Nalpha"])             <- "numeric"
	mode(solution.df[,"Sn,alpha"])               <- "numeric"
	mode(solution.df[,"z"])                      <- "numeric"
	mode(solution.df[,"(z-Sn,alpha)"])           <- "numeric"
	mode(solution.df[,"((z-Sn,alpha)+)"])        <- "numeric"
	mode(solution.df[,"((z-Sn,alpha)^(n-1))"])   <- "numeric"
	mode(solution.df[,"((-1)^Nalpha)*((z-Sn,alpha)^(n-1))"]) <- "numeric"
	solution.df[,"(z-Sn,alpha)"]           <- as.numeric(solution.mat[,"z"])-as.numeric(solution.mat[,"Sn,alpha"])
	solution.df[,"((z-Sn,alpha)+)"]        <- (solution.df[,"(z-Sn,alpha)"] + abs(solution.df[,"(z-Sn,alpha)"]))/2
	if(cumulative){
		solution.df[,"((z-Sn,alpha)^(n-1))"]               <- solution.df[,"((z-Sn,alpha)+)"]^(n)
	} else {
		solution.df[,"((z-Sn,alpha)^(n-1))"]               <- solution.df[,"((z-Sn,alpha)+)"]^(n-1)
	}
	solution.df[,"((-1)^Nalpha)*((z-Sn,alpha)^(n-1))"] <- as.numeric(solution.mat[,"(-1)Nalpha"])*solution.df[,"((z-Sn,alpha)^(n-1))"]
	part2  <- sum(solution.df[,"((-1)^Nalpha)*((z-Sn,alpha)^(n-1))"])
	result <- part1*part2
	result
}

#' @title PDF or CDF for the product of n independent uniformly continuous random variables.
#' 
#' Evaluates the probability density function or cumulative distribution of Z=z for Z(z) = Xi*...*Xn, for i = 1...n, Xi ~ U(ai,bi), where 0 < ai < bi
#' The PDF and CDF is calculated using Equations 3.4 and 3.9, respectively, from Ishihara, Tatsuo. 2002. The Distribution of the Sum and the Product of Independent Uniform Random Variables Distributed at Different Intervals. Transactions - Japan Society for Industrial and Applied Mathematics, 12(3), 197-208.
#' 
#' @param z Number or numerical vector of quantiles
#' @param a Vector holding the mins of each X
#' @param b Vector holding the max of each X
#' @param cumulative Whether or not return the cumulative probability instead of the probability. Default is FALSE. This parameter does not yet calculate CDF correctly.
#' @return Number or numerical vector with either the probability density or the cumulative density for each input value z.
#' @export produnif
produnif <- function(z,a,b,cumulative=F){
	if(!cumulative & (z <= prod(a) | z > prod(b))){
		return(0)
	}
	if(cumulative & z <= prod(a)){
		return(0)
	}
	if(cumulative & z >= prod(b)){
		return(1)
	}
	C  <- b-a
	CN <- prod(C)
	n  <- length(a)
	powerSet.list                    <- pset(1:n)
	solution.mat                     <- matrix(data=NA,nrow=length(powerSet.list),ncol=5)
	colnames(solution.mat)           <- c("z","No.","alpha","(-1)Nalpha","Pn,alpha")
	solution.df                      <- as.data.frame(solution.mat)
	mode(solution.df[,"z"])          <- "numeric"
	mode(solution.df[,"No."])        <- "factor"
	mode(solution.df[,"(-1)Nalpha"]) <- "numeric"
	mode(solution.df[,"Pn,alpha"])   <- "numeric"
	solution.df[,"z"]                <- rep(z,nrow(solution.df))
	solution.df[,"No."]              <- 1:nrow(solution.df)
	solution.df[,"alpha"]            <- sapply(X = powerSet.list, FUN=paste, collapse = " ")
	solution.df[,"(-1)Nalpha"]       <- (-1)^lengths(powerSet.list)
	if(!cumulative){
		part1 <- 1/((factorial(n-1))*CN)
		for(i in 1:nrow(solution.df)){
			alpha.i <- rep(0,n)
			if(i>1){
				alpha.i[1:n %in% powerSet.list[[i]]] <- 1
			}
			Pn.i    <- prod(a + (alpha.i*C))
			solution.df[i,"Pn,alpha"] <- Pn.i
		}
		solution.df[,"ln(x/Pn,alpha)+"]      <- (log(z/solution.df[,"Pn,alpha"]) + abs(log(z/solution.df[,"Pn,alpha"])))/2
		solution.df[,"(ln(x/Pn,alpha)^(d))"] <- solution.df[,"ln(x/Pn,alpha)+"]^(n-1)
		result <- (1/((factorial(n-1))*CN))*sum(solution.df[,"(-1)Nalpha"]*solution.df[,"(ln(x/Pn,alpha)^(d))"])
	} else {
		for(i in 1:nrow(solution.df)){
			alpha.i <- rep(0,n)
			if(i>1){
				alpha.i[1:n %in% powerSet.list[[i]]] <- 1
			}
			Pn.i    <- prod(a + (alpha.i*C))
			solution.df[i,"Pn,alpha"] <- Pn.i
			part2.i <- solution.df[i,"(-1)Nalpha"]
			part3.j <- list(); length(part3.j) <- length(1:(n-1))
			for(j in 1:(n-1)){
				part3.j1 <- ((-1)^(j+1))/(factorial(n-j))
				part3.j2 <- ((log((z/Pn.i)) + abs(log((z/Pn.i))))/2)^(n-j)
				#part3.j[[j]] <- (part3.j1*part3.j2)+(part3.j3*part3.j4)
				#part3.j[[j]] <- (part3.j1*part3.j2)+((-1)^(n-1))*(((z-Pn.i)+abs(z-Pn.i))/2)
				part3.j[[j]] <- (part3.j1*part3.j2)
			}
			part3.i <- z*sum(unlist(part3.j))
			solution.df[i,"part3.i"] <- part3.i
			part4.i <- ((-1)^(n-1))*(((z-Pn.i)+abs(z-Pn.i))/2)
			solution.df[i,"part4.i"] <- part4.i
		}
		solution.df[,"parts.i"] <- solution.df[,"(-1)Nalpha"]*(solution.df[,"part3.i"]+solution.df[,"part4.i"])
		result <- (sum(solution.df[,"parts.i"])/CN)
	}
	result
}

#' @title plot PDF for sum or product of two uniform distributions
#' 
#' Plots a line chart for the probability density distribution (PDF) or cumulative distribution (CDF) for the sum or product of uniform distributions on positive Reals.
#' Can handle an arbitrary number of input distributions, which must be independent, and can have different ranges.
#'
#' @param a A numerical vector holding the minimum values of input uniform distributions.
#' @param b A numerical vector holding the maximum values of input uniform distributions.
#' @param FUN The function to use to determine PDF or CDF; one of 'simple', 'sum', or 'product'. If a and b are length 1 vectors, 'simple' is used.
#' @param cumulative If TRUE, the CDF will be plotted instead of the PDF
#' @param buffer fraction of extra, zero-probability region along x axis to plot
#' @param show.xy.pd When FUN = 'sum' or 'product', if 'show.xy.pd' = TRUE the pdf or cdf of component distributions will be plotted along with their sum or product.
#' @param line.lty type of curve
#' @param line.col color of curve
#' @param add Whether or not to the graph should be added to the existing plot. Default is FALSE.
#' @param xlim min and max limits for the x-axis plotting area
#' @param ylim min and max limits for the y-axis plotting area
#' @return A line graph showing the PDF of the sum or product of the two input uniform distributions
#' @export plot.funUnif
plot.funUnif <- function(a,b,FUN="sum",cumulative=F,buffer=0.1,nZ=500,show.xy.pd=F,line.col=c("red"),line.lty=c(1),add=F,xlim=NULL,ylim=NULL){
	if(FUN=="simple" | (length(a) == 1 & length(b) == 1)){
		#a=4e+05;b=9e+05;buffer=0.1;
		x0=(a-abs((a*buffer)))
		x1=(b+abs((b*buffer)))
		if(is.null(xlim)){
			xlim=c(x0,x1)
		}
		z.samples  <- seq(from=x0,to=x1,length.out=nZ)
		z.pd       <- sapply(X=z.samples,FUN=function(x){dunif(x=x,min=a,max=b)})
		if(is.null(ylim)){
			ylim=c(0,(max(z.pd)*1.05))
		}
		if(!add){
			plot(z.samples,z.pd,col="white",xlab="",ylab="",main="",axes=F,xlim=xlim)
			axis(side=1,at=xlim,labels=T)
		}
		lines(z.samples,z.pd,col=line.col[1],lty=line.lty[1])
	} else {
		if(FUN=="sum"){
			zmin = sum(a)
			zmax = sum(b)
		}
		if(FUN=="product"){
			zmin = prod(a)
			zmax = prod(b)
		}
		if(!FUN %in% c("sum","product")){
			stop("'FUN' must be 'sum' or 'product'")
		}
		zmin.buffered <- zmin - abs((zmax-zmin)*buffer)
		zmax.buffered <- zmax + abs((zmax-zmin)*buffer)
		a2 <- a - abs((b-a)*buffer)
		b2 <- b + abs((b-a)*buffer)
		zrange <- c(zmin.buffered,zmax.buffered)
		xrange <- c(a2[1],b2[1])
		yrange <- c(a2[2],b2[2])
		z.samples <- seq(from=zrange[1],to=zrange[2],length.out=nZ)
		if(FUN=="sum"){
			z.pd <- sapply(z.samples,FUN=sumunif,a=a,b=b,cumulative=cumulative)
		}
		if(FUN=="product"){
			z.pd <- sapply(z.samples,FUN=produnif,a=a,b=b,cumulative=cumulative)
		}
		if(is.null(ylim)){
			ylim <- c(0,(max(z.pd)*1.05))
		}
		if(show.xy.pd){
			if(length(a)>2){
				stop("'Show.xy.pdf' must be FALSE when using more than two input distributions")
			}
			x.samples <- seq(from=a2[1],to=b2[1],length.out=nZ)
			y.samples <- seq(from=a2[2],to=b2[2],length.out=nZ)
			x.pd      <- dunif(x.samples,min=a[1],max=b[1])
			y.pd      <- dunif(y.samples,min=a[2],max=b[2])
			if(is.null(xlim)){
				xlim <- c(min(a2,zrange[1]), max(b2,zrange[2]))
			}
			if(!add){
				plot(c(x.samples,y.samples,z.samples),c(x.pd,y.pd,z.pd),col="white",xlab="",ylab="",main="",axes=F,xlim=xlim,ylim=ylim)
				axis(side=1,at= xlim, labels=T)
			}
			if(length(line.col) < 3){
				line.col=c("red","blue","green")
			}
			if(length(line.lty) < 3){
				line.lty=c(1,2,2)
			}
			lines(z.samples,z.pd,col=line.col[1],lty=line.lty[1])
			lines(x.samples,x.pd,col=line.col[2],lty=line.lty[2])
			lines(y.samples,y.pd,col=line.col[3],lty=line.lty[3])
		} else {
			if(!add){
				if(is.null(xlim)){
					xlim <- zrange
				}
				plot(z.samples,z.pd,col="white",xlab="",ylab="",main="",axes=F,xlim=xlim,ylim=ylim)
				axis(side=1,at=xlim,labels=T)
			}
			lines(z.samples,z.pd,col=line.col[1],lty=line.lty[1])
		}
	}
}

#' @title PDF sum of normal distributions
#' 
#' Probability distribution for the sum of n normal distributions
#' 
#' @param x sample quantile
#' @param m numerical vector with means of the normal distributions
#' @param v numerical vector with variances of the normal distributions
#' @return probability distribution at x for the sum of the normal distributions with means m and variances v
#' @export sumnorm.PDF
sumnorm.PDF <- function(x,m,v){
	dist <- distr::Norm(sum(m),sqrt(sum(v)))
	return(distr::d(dist)(x))
}



##### Function to return the powerset for a numeric or character vector of set elements.
##' 
#pset <- function (x) {
#	m   <- length(x)
#	out <- list(x[c()])
#	for (i in seq_along(x)) {
#		out <- c(out, lapply(X=out[lengths(out) < m], FUN=function(y){c(y,x[i])}))
#	}
#	out
#}




