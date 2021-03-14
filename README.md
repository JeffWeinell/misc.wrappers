# misc.wrappers
 Collection of wrappers for programs that I did not create.

```
library(BiocManager)
BiocManager::install("JeffWeinell/misc.wrappers",auth_token="323d9e4cd00247a39a805dbb66f37db6403cfb8b")
```



<!--
```
library(ade4)
library(ape)
library(stringr)
library(vcfR)
library(adegenet)
library(pegas)
# library(poppr)

# BiocManager::install(c("dplyr"))
# BiocManager::install(c("spData"))
# BiocManager::install(c("sf","spdep","adegenet"))
# 
# BiocManager::install(c("quadprog","phangorn"))
# 
# quadprog
# 
### Need to do module load gdal before installing sf. That didnt work though.
# module load compiler/gcc/8.3
# module load gdal/3.0.0
# module load geos/3.7.2
# module load proj/6.0.0
# module load R/4.0
# R
# .libPaths("/panfs/pfs.local/home/j926w878/programs/R-packages")
# old_path <- Sys.getenv("PATH")
# Sys.setenv(PATH = paste(old_path, "/panfs/pfs.local/software/7/install/modulefiles/Core", sep = ":"))
# install.packages("sf")
# 
# #BiocManager::install("ape")
# #BiocManager::install("ade4")
# install.packages("pegas")
# #install.packages("stringr")
# install.packages("adegenet")
# 
### eems for snp data requires three files (.diffs, .coord, and .outer). This liitle tutorial will help you make the first 2/3 of those files and alter the third.
## this first part will help you make the .diffs file
# read in 1 snp per locus MAF
# vmaf        <- vcfR::read.vcfR("/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/Oxyrhabdium_populations.snps_28Feb2021.vcf")
# ## Drop the 63rd individual, "Oxyrhabdium-modestum-widespread_KU345787_Siniloan-Laguna_15", because it has almost no data.
# vmaf        <- vmaf[,,c(1:62,64:65)]
# 
# vmaf_Oxyrhabdium <- vmaf
# vmaf_leporinum   <- vmaf[,,1:29]
# vmaf_leporinum_Luzon <- vmaf_leporinum[,,c(1:3,6:11,13:29)]
# vmaf_modestum    <- vmaf[,,30:42]
# vmaf_cf.modestum <- vmaf[,,c(43:65)]
# vmaf_cf.modestum <- vmaf_cf.modestum[,,c(1:20,22:23)]
#
### Convert to genind object
# gen_leporinum_Luzon <- vcfR::vcfR2genind(vmaf_leporinum_Luzon)
# ploidy(gen) <- 2
# 
# ### From https://github.com/dipetkov/eems/blob/master/str2diffs/README.md
# #data       <- gen
# stopifnot(identical(gen@type, 'codom'))
# Geno       <- gen@tab
# ## Names of loci that are not biallelic
# multi.loci <- names(which(data@loc.n.all != 2))
# ## Explanation: 
# ## Suppose we want to remove locus, `L100` with alleles, `L100.00`, `L100.01`, `L100.02`, then detect columns whose names matches the regular expression `^L100\\.\\d+$`
# ## Columns corresponding to the non-biallelic loci.
# # multi.cols <- which(grepl(paste0("^", multi.loci, "\\.\\d+$", collapse = "|"), colnames(Geno)))
# ### Faster version
# multi.cols <- gsub("\\..+","",colnames(Geno)) %in% multi.loci
# ### Removing columns for loci that are not biallelic
# if(length(which(multi.cols))){
#   Geno <- Geno[, -which(multi.cols)]
# }
# dim(Geno)
### Two ways to do bed2diffs
# bed2diffs_v2 <- function(Geno) {
#     nIndiv <- nrow(Geno)
#     nSites <- ncol(Geno)
#     Miss   <- is.na(Geno)
#     ## Impute NAs with the column means (= twice the allele frequencies)
#     Mean   <- matrix(colMeans(Geno, na.rm = TRUE), ## a row of means
#                    nrow = nIndiv, ncol = nSites, byrow = TRUE) ## a matrix with nIndiv identical rows of means
#     ## Set the means that correspond to observed genotypes to 0
#     Mean[Miss == 0] <- 0
#     ## Set the missing genotypes to 0 (used to be NA) 
#     Geno[Miss == 1] <- 0
#     Geno     <- Geno + Mean
#     ## Compute similarities
#     Sim      <- Geno %*% t(Geno) / nSites
#     ## self-similarities
#     SelfSim  <- diag(Sim)
#     ## vector of 1s
#     vector1s <- rep(1, nIndiv)
#     ## This chunk generates a `diffs` matrix
#     Diffs    <- SelfSim %*% t(vector1s) + vector1s %*% t(SelfSim) - 2 * Sim
#     Diffs
# }
# 
# bed2diffs_v1 <- function(Geno) {
#   nIndiv <- nrow(Geno)
#   nSites <- ncol(Geno)
#   Diffs <- matrix(0, nIndiv, nIndiv)
#   
#   for (i in seq(nIndiv - 1)) {
#     for (j in seq(i + 1, nIndiv)) {
#       x <- Geno[i, ]
#       y <- Geno[j, ]
#       Diffs[i, j] <- mean((x - y)^2, na.rm = TRUE)
#       Diffs[j, i] <- Diffs[i, j]
#     }
#   }
#   Diffs
# }
```
-->

```
# ### Create diffs matrices
# Diffs_v1 <- round(bed2diffs_v1(Geno), digits = 6)
# Diffs_v2 <- round(bed2diffs_v2(Geno), digits = 6)
# 
# ### Pick the diff file that has only 1 positive eigenvalue. JLW: What if neither have only one??? There should be as many eigenvalues as the length of the diagonal of the input matrix.
# eigenvalues.Diffs_v1 <- sort(round(eigen(Diffs_v1)$values, digits = 2))
# eigenvalues.Diffs_v2 <- sort(round(eigen(Diffs_v2)$values, digits = 2))
# 
# # write.table(Diffs_v1, "/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/leporinum_Diffs1.diffs", col.names = FALSE, row.names = FALSE, quote = FALSE)
# # write.table(Diffs_v2, "/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/leporinum_Diffs2.diffs", col.names = FALSE, row.names = FALSE, quote = FALSE)

## Create diffs matrices.
Diffs_leporinum_Luzon <- genind2diffs(genind.obj=gen_leporinum_Luzon,output.file="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/leporinum_Luzon.diffs")

#write.table(Diffs_v1, "/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/leporinum-Luzon_Diffs1.diffs", col.names = FALSE, row.names = FALSE, quote = FALSE)
#write.table(Diffs_v2, "/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/leporinum-Luzon_Diffs2.diffs", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Label each individual as a pop so that you can use dist.genpop to estimate distances properly (this is for 238 individuals)
# pop     <- as.factor(seq(1,238,1))
# pop     <- pop(gen)
# gen_pop <- genind2genpop(gen)
# gendist <- dist.genpop(gen_pop, method=5, upper=T, diag=T)
# write.table(as.matrix(gendist),"rats.diffs", col.names=F,row.names=F)

### Make the .coord file. Make sure it and the .demes file are the in the same orientation (e.g., long lat). the .demes file is made with the following link: 
### check names and lats and lons

#lats  <- read.table("/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/Oxyrhabdium_snps_LEA_v2_coordinates.txt",sep="\t")[,c(3,2)]
#colnames(lats) <- c("long","lat")
#z     <- lats

# lats  <- read.csv( "East_FL_lats.csv",header=T)
# names <- gen$ind.names
# names <- rownames(Geno)
# geo   <- cbind(names,lats)
# z     <- read.csv( "East_FL_lats.csv",header=T)

#make sure long/lat are first to match up to the .outer outline file
# z2    <- cbind(z$long,z$lat)

# write.table(z2, "/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/Oxyrhabdium.coord", row.names=F,col.names=F)
# Oxyrhadbium.coord     <- read.table("/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/Oxyrhabdium.coord", header=F)
# leporinum.coord       <- Oxyrhadbium.coord[1:29,]
# leporinum_Luzon.coord <- leporinum.coord[c(1:3,6:11,13:29),]
# write.table(leporinum.coord, "/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/leporinum.coord", row.names=F,col.names=F,sep=" ",quote=F)
# write.table(leporinum_Luzon.coord, "/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/leporinum_Luzon.coord", row.names=F,col.names=F,sep=" ",quote=F)
```




```
.libPaths("/panfs/pfs.local/home/j926w878/programs/R-packages")
library(vcfR)

vmaf                 <- vcfR::read.vcfR("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/Oxyrhabdium_populations.snps_28Feb2021.vcf")
vmaf                 <- vmaf[,,c(1:62,64:65)]
vmaf_Oxyrhabdium     <- vmaf
vmaf_leporinum       <- vmaf[,,1:29]
vmaf_leporinum_Luzon <- vmaf_leporinum[,,c(1:3,6:11,13:29)]
vmaf_modestum        <- vmaf[,,30:42]
vmaf_cf.modestum     <- vmaf[,,43:64]
gen_leporinum_Luzon  <- vcfR::vcfR2genind(vmaf_leporinum_Luzon)

test <- runEEMs(eems.path = "/panfs/pfs.local/home/j926w878/programs/eems-master/runeems_snps/src/runeems_snps", input.data="/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/leporinum_Luzon.diffs",path.coord="/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/leporinum_Luzon.coord",path.outer="/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/leporinum_Luzon.outer",output.path="/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_leporinum_Luzon",setup.only=T,n.sites=27436)

#
eems.path = "/panfs/pfs.local/home/j926w878/programs/eems-master/runeems_snps/src/runeems_snps"

#################
### Ahaetulla ###
Ahaetulla.vcf <- vcfR::read.vcfR("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/Ahaetulla_snps.vcf")
## Remove individual with mostly missing data, and individual from Samar (Because only looking at Luzon right now).
Ahaetulla_Luzon.vcf   <- Ahaetulla.vcf[,,c(1:7,9:16,18:21)]
Ahaetulla_Luzon.gen   <- vcfR::vcfR2genind(Ahaetulla_Luzon.vcf)
Ahaetulla_Luzon.diffs <- genind2diffs(genind.obj=Ahaetulla_Luzon.gen,output.file="/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/Ahaetulla_data.diffs")
# Move the diffs, coords, and outer file to the cluster and do the next line on the cluser.
Ahaetulla.eems  <- runEEMs(eems.path = "/panfs/pfs.local/home/j926w878/programs/eems-master/runeems_snps/src/runeems_snps", input.data="/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/Ahaetulla_data.diffs",path.coord="/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/Ahaetulla_v2_Luzon.coord",path.outer="/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/old_EEMS_Calamaria-gervaisii_Luzon/data/data.outer",output.path="/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_Ahaetulla-prasina_Luzon",setup.only=T,n.sites=19670)

#################
### Calamaria ###
# First work locally to generate diffs file.
Calamaria.vcf   <- vcfR::read.vcfR("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/Calamaria_snps.vcf")
Calamaria.gen   <- vcfR::vcfR2genind(Calamaria.vcf)
Calamaria.diffs <- genind2diffs(genind.obj=Calamaria.gen,output.file="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/EEMS/EEMS_Calamaria-gervaisii_Luzon/data/Calamaria_data.diffs")
# Move the diffs, coords, and outer file to the cluster and do the next line on the cluser.
Calamaria.eems  <- runEEMs(eems.path = "/panfs/pfs.local/home/j926w878/programs/eems-master/runeems_snps/src/runeems_snps", input.data="/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/Calamaria_data.diffs",path.coord="/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/old_EEMS_Calamaria-gervaisii_Luzon/data/data.coord",path.outer="/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/old_EEMS_Calamaria-gervaisii_Luzon/data/data.outer",output.path="/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_Calamaria-gervaisii_Luzon",setup.only=T,n.sites=15166)

#gen_leporinum_Luzon <- vcfR::vcfR2genind(vmaf_leporinum_Luzon)
#/panfs/pfs.local/home/j926w878/programs/eems-master/runeems_snps/src/runeems_snps --params /panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_leporinum_Luzon/params/params-chain1_300.ini
#/panfs/pfs.local/home/j926w878/programs/eems-master/runeems_snps/src/runeems_snps --params /panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_leporinum_Luzon/params/params-chain2_300.ini
#/panfs/pfs.local/home/j926w878/programs/eems-master/runeems_snps/src/runeems_snps --params /panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_leporinum_Luzon/params/params-chain3_300.ini

########
### Plotting results. First using an example.
library(reemsplots2)
library(ggplot2)
library(gridExtra)
theme(plot.margin=grid::unit(c(1,1,1,1), "inches"))
plots.chain1 <- suppressWarningsmake_eems_plots("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_leporinum_Luzon2/mcmc/300Demes-chain1", longlat = TRUE))
plots.chain2 <- suppressWarnings(make_eems_plots("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_leporinum_Luzon2/mcmc/300Demes-chain2", longlat = TRUE))
plots.chain3 <- suppressWarningsmake_eems_plots("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_leporinum_Luzon2/mcmc/300Demes-chain3", longlat = TRUE))
ggsave("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_leporinum_Luzon2/plots_chain1.pdf", gridExtra::marrangeGrob(grobs = plots.chain1, nrow=1, ncol=1),width=11,height=8.5)
ggsave("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_leporinum_Luzon2/plots_chain2.pdf", gridExtra::marrangeGrob(grobs = plots.chain2, nrow=1, ncol=1),width=11,height=8.5)
ggsave("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_leporinum_Luzon2/plots_chain3.pdf", gridExtra::marrangeGrob(grobs = plots.chain3, nrow=1, ncol=1),width=11,height=8.5)

# pdf(file="/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_leporinum_Luzon2/plots_chain1.pdf",width=11,height=8)
# Using mfrow doesnt work because these plots are class ggplot
# par(mfrow=c(1,2))
# plot(plots[[1]])
# plot(plots[[2]])
# plot(plots[[3]])
# plot(plots[[4]])

# grid.arrange(plots[[1]], plots[[2]], ncol=2)
# ggsave("/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/example_output.pdf", arrangeGrob(plots[[1]], plots[[2]]))
### This seems to work but adds an extra blank page at the start.
# ggsave("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_leporinum_Luzon2/plots_chain1.pdf", gridExtra::marrangeGrob(grobs = plots, nrow=1, ncol=2))


.libPaths("/panfs/pfs.local/home/j926w878/programs/R-packages")
library(reemsplots2)
library(ggplot2)
library(gridExtra)
### Generate EEMS plots
plots.chain1 <- suppressWarnings(make_eems_plots("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_Calamaria-gervaisii_Luzon/mcmc/300Demes-chain1", longlat = TRUE))
plots.chain2 <- suppressWarnings(make_eems_plots("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_Calamaria-gervaisii_Luzon/mcmc/300Demes-chain2", longlat = TRUE))
plots.chain3 <- suppressWarnings(make_eems_plots("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_Calamaria-gervaisii_Luzon/mcmc/300Demes-chain3", longlat = TRUE))
### Update theme margins attribute
plots.chain1.v2 <- lapply(X=plots.chain1,FUN=function(x){x+theme(plot.margin = unit(c(0.25,0.25,0.25,0), "inches"))})
plots.chain2.v2 <- lapply(X=plots.chain2,FUN=function(x){x+theme(plot.margin = unit(c(0,0.25,0.25,0), "inches"))})
plots.chain3.v2 <- lapply(X=plots.chain3,FUN=function(x){x+theme(plot.margin = unit(c(0,0.25,0.25,0.25), "inches"))})
### merge the lists of ggplots into a single list of ggplots with ith plot of each chain next to each other
# plots.all <- do.call(c,lapply(X=c(1:8),FUN=function(x){list(plots.chain1.v2[[x]],plots.chain2.v2[[x]],plots.chain3.v2[[x]])}))
plots.all <- c(plots.chain1.v2,plots.chain2.v2,plots.chain3.v2)
###
# grob.all  <- gridExtra::marrangeGrob(grobs = plots.all, nrow=1, ncol=3)
grob.all  <- gridExtra::arrangeGrob(grobs = plots.all, nrow=8, ncol=3)
###
res <- pdf("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_Calamaria-gervaisii_Luzon/plots_all-chains.pdf",width=20,height=40)
lapply(X=grob.all,FUN=plot)
dev.off()

### Different way of storing the plots so that the pdf function can be used
# Important to run these lines before running the pdf line, otherwise a blank page will appear at the start.
grob.chain1 <- gridExtra::marrangeGrob(grobs = plots.chain1.v2, nrow=1, ncol=1)
grob.chain2 <- gridExtra::marrangeGrob(grobs = plots.chain2.v2, nrow=1, ncol=1)
grob.chain3 <- gridExtra::marrangeGrob(grobs = plots.chain3.v2, nrow=1, ncol=1)
# Now run pdf, lapply the plots, and then dev.off
pdf("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_Calamaria-gervaisii_Luzon/plots_chain1.pdf",width=11,height=8.5)
lapply(X=grob.chain1,FUN=plot)
dev.off()
pdf("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_Calamaria-gervaisii_Luzon/plots_chain3.pdf",width=11,height=8.5)
lapply(X=grob.chain2,FUN=plot)
dev.off()
pdf("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_Calamaria-gervaisii_Luzon/plots_chain3.pdf",width=11,height=8.5)
lapply(X=grob.chain3,FUN=plot)
dev.off()


#######
#!/bin/bash
# sbatch --nodes=1 --ntasks-per-node=30 --time=6:00:00 --partition=sixhour /panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_leporinum_Luzon2/ploteems2.sh
.libPaths("/panfs/pfs.local/home/j926w878/programs/R-packages")
library(reemsplots2)
library(ggplot2)
library(gridExtra)
### Load a more detailed outline of Luzon than the polygon used for running eems
Luzon.outer           <- read.table("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/Luzon.outer",sep=" ")
colnames(Luzon.outer) <- c("x","y")
### Paths to mcmc output
mcmcpath.chain1 <- "/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_Calamaria-gervaisii_Luzon/mcmc/300Demes-chain1"
mcmcpath.chain2 <- "/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_Calamaria-gervaisii_Luzon/mcmc/300Demes-chain2"
mcmcpath.chain3 <- "/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_Calamaria-gervaisii_Luzon/mcmc/300Demes-chain3"
### Generate EEMS plots
plots.chain1 <- suppressWarnings(make_eems_plots(mcmcpath.chain1, longlat = TRUE))
plots.chain2 <- suppressWarnings(make_eems_plots(mcmcpath.chain2, longlat = TRUE))
plots.chain3 <- suppressWarnings(make_eems_plots(mcmcpath.chain3, longlat = TRUE))
### Update theme margins attribute
plots.chain1.v2 <- lapply(X=plots.chain1,FUN=function(x){x+theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "inches"))})
plots.chain2.v2 <- lapply(X=plots.chain2,FUN=function(x){x+theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "inches"))})
plots.chain3.v2 <- lapply(X=plots.chain3,FUN=function(x){x+theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "inches"))})
### Add Luzon.outer geometry path for the first 4 ggplot object
plots.chain1.v2[1:4] <- lapply(X=plots.chain1.v2[1:4],FUN=function(x){x+theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "inches"))+geom_path(data=Luzon.outer)})
plots.chain2.v2[1:4] <- lapply(X=plots.chain2.v2[1:4],FUN=function(x){x+theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "inches"))+geom_path(data=Luzon.outer)})
plots.chain3.v2[1:4] <- lapply(X=plots.chain3.v2[1:4],FUN=function(x){x+theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "inches"))+geom_path(data=Luzon.outer)})

### merge the lists of ggplots into a single list of ggplots with ith plot of each chain next to each other
plots.all <- do.call(c,lapply(X=c(1:8),FUN=function(x){list(plots.chain1.v2[[x]],plots.chain2.v2[[x]],plots.chain3.v2[[x]])}))
###
# grob.all  <- gridExtra::marrangeGrob(grobs = plots.all, nrow=1, ncol=3)
grob.all  <- gridExtra::arrangeGrob(grobs = plots.all, nrow=8, ncol=3)
pdf("/panfs/pfs.local/home/j926w878/work/ddRAD/EEMS/EEMS_Calamaria-gervaisii_Luzon/plots_all-chains.pdf",width=21,height=42)
plot(grob.all)
dev.off()

####### Attempting to plot an outline for Luzon on top of the EEMS output
library(sp)
Luzon.coords   <- read.table("/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/leporinum_Luzon.outer",sep=" ")
Luzon.polygon  <- Polygon(Luzon.coords,hole=F)
Luzon.polygons <- Polygons(list(Luzon.polygon),ID="Luzon")
Luzon.sp       <- SpatialPolygons(Srl=list(Luzon.polygons))
plot(Luzon.sp)

outer.coords    <- read.table("/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/leporinum_Luzon2.outer",sep=" ")
outer.polygon  <- Polygon(outer.coords,hole=F)
outer.polygons <- Polygons(list(outer.polygon),ID="Luzon.outer")
outer.sp       <- SpatialPolygons(Srl=list(outer.polygons))

### this is the first ggplot for chain2 Ahaetulla
plot(plots.chain2[[1]])

###### This works to show Luzon outline!!!
Luzon.coords           <- read.table("/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/leporinum_Luzon.outer",sep=" ")
colnames(Luzon.coords) <- c("x","y")
plots.chain2[[1]] + geom_path(data=Luzon.coords)
plots.chain2.v2 <- lapply(X=plots.chain2,FUN=function(x){x+theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "inches"))+geom_path(data=Luzon.outer})

# required packages
# library(sp)                 ### Used for everything
# library(maps)               ### map.axes function
# library(rgeos)              ### gBuffer and gArea function
# library(rgdal)              ### showWKT function
# library(raster)             ### crs and bind function
# library(smoothr)            ### densify function for adding outer points by interpolating along perimeter
# library(sampSurf)           ### spCircle function
# library(geosphere)          ### perimeter and areaPolygon functions
# library(alphahull)          ### ahull function
# library(adehabitatHR)       ### function to generate minimum convex polygon
# library(rnaturalearth)      ### Global geographic vectors
# library(rnaturalearthhires) ### high resolution data for rnaturalearth package


#### Include metrics describing polygon shape:
## Described here: http://www.umass.edu/landeco/teaching/landscape_ecology/schedule/chapter9_metrics.pdf
# PARA = perimeter-area ratio; confounded with polygon size, such that a larger polygon with the same shape will have a smaller PARA.
# SHAPE = shape index; normalized PARA, in which the complexity of patch shape is compared to a standard shape (square) of the same size, thereby alleviating the size dependency problem of PARA. Values > 1 represent increasing departure from square shape.
# FRAC = fractal dimension index; a normalized shape index based on PARA, in which the perimeter and area are log transformed; FRAC = D = (2*log(P))/(log(A)), where P = perimeter and A = area. For simple Euclidean shapes (e.g., circles, squares), D = 1, and D approaches 2 as the polygon becomes increassingly complex.
# PAFRAC = perimeter-area fractal dimension; a similar index to FRAC, but applied to a collection of patches at the class or landscape level.
# CIRCLE = circumscribing circle index; based on the ratio of patch area to the area of the smallest circumscribing circle; providing a measure of overall patch elongation. A highly convoluted but narrow patch will have a relatively low related circumscribing circle index due to the relative compactness of the patch. Conversely, a narrow and elongated patch will have a relatively high circumscribing circle index. This index may be particularly useful for distinguishing patches that are both linear (narrow) and elongated.
# CONTIG = contiguity index; method of assessing patch shape based on the spatial connectedness or contiguity of cells within a patch; large contiguous patches will result in larger contiguity index values.

```



<!--
test <- create.outer(coords="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/Oxyrhabdium.coord")
test <- create.outer(coords="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/Ahaetulla_v2.coord")
test <- create.outer(coords="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/Ahaetulla_v2_Luzon.coord")
test <- create.outer(coords="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/leporinum_Luzon.coord")
test <- create.outer(coords="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/leporinum.coord",max.fractal.dimension=1.5)
test <- create.outer(coords="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/EEMS_Katie/mainland_maculilabris.coord",method=1,coords.radius=0.5)
#
#
coords="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/Oxyrhabdium.coord";buffer.adj=0;coords.radius=0.01;method=1;plot.outer=T;counter.clockwise=T;output.path=NULL
coords="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/Ahaetulla_v2.coord";buffer.adj=0;coords.radius=0.01;method=1;plot.outer=T;counter.clockwise=T;output.path=NULL
coords="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/Ahaetulla_v2_Luzon.coord";buffer.adj=0;coords.radius=0.01;method=1;plot.outer=T;counter.clockwise=T;output.path=NULL
coords="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/leporinum_Luzon.coord";buffer.adj=0;coords.radius=0.01;method=1;plot.outer=T;counter.clockwise=T;output.path=NULL
coords="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/leporinum.coord";buffer.adj=0;coords.radius=0.01;method=3;plot.outer=T;counter.clockwise=T;output.path=NULL;max.fractal.dimension=1.5
coords="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/EEMS_Katie/mainland_maculilabris.coord";buffer.adj=0;coords.radius=2;method=1;plot.outer=T;counter.clockwise=T;output.path=NULL
coords="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/EEMS/EEMS_Katie/mainland_maculilabris.coord";buffer.adj=0;method=2;plot.outer=T;counter.clockwise=T;output.path=NULL
#
#
spdf_Brazil_10 <- rnaturalearth::ne_countries(scale=10,country="Brazil")
Brazil.points.sampled <- spsample(spdf_Brazil_10, n=30, type="random")
Brazil.coords.sampled <- sp2coords(Brazil.points.sampled)
Brazil.test.random    <- create.outer(Brazil.coords.sampled)
-->
