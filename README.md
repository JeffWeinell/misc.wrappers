# misc.wrappers
 **Under Developement**

### Package Description
This package includes functions to fascilitate running some popular programs used for assessing population structure, demographic modelling, and genetic/genomic and phylogenetic/genomic programs, and to summarize results in useful ways. The main functions are shown in the table below along with example input and output files. The values used for function arguments, not shown in table, also control analyses and should be set mindfully - see examples (will add soon) of function usage.

function | description | example input files | example output files
---|---|---|-----
```plot_model``` | Makes a cartoon figure of a model described in a pair of .tpl and .est files | [model10_K3.tpl](inst/extdata/example.tpl), [model10_K3.est](inst/extdata/example.est) | [model10_K3.pdf](inst/extdata/example_model.pdf)
```create.outer``` | Generates a set of outer habitat coordinates that can be used for EEMS | [coords.txt](inst/extdata/createouter_exampleInput_coords.txt) | method 1: [outer1.txt](inst/extdata/createouter_exampleOutput_method1_outer.text), [outer1.pdf](inst/extdata/createouter_exampleOutput_method1_outer.pdf)<br /><br /> method 2: [outer2.txt](inst/extdata/createouter_exampleOutput_method2_outer.text), [outer2.pdf](inst/extdata/createouter_exampleOutput_method2_outer.pdf)<br /><br /> method 3: [outer3.txt](inst/extdata/createouter_exampleOutput_method3_outer.text), [outer3.pdf](inst/extdata/createouter_exampleOutput_method3_outer.pdf)
```sim.vcf ``` | Simulate SNPs and save as a VCF; optionally introduce missing data; optionally simulate locality data | [example.vcf.gz](inst/extdata/example.vcf.gz) (optional) | K=3, allopatric: [simK3.vcf.gz](inst/extdata/simK3.vcf.gz), [simK3_coords.txt](inst/extdata/simK3_coords.txt), [simK3_coords_map.pdf](inst/extdata/simK3_coords_map.pdf)<br /><br /> K=4, some contact: [simK4.vcf.gz](inst/extdata/simK4.vcf.gz), [simK4_coords.txt](inst/extdata/simK4_coords.txt), [simK4_coords_map.pdf](inst/extdata/simK4_coords_map.pdf)
```run_DAPC``` | Pipeline for running DAPC and graphing results | [simK4.vcf.gz](inst/extdata/simK4.vcf.gz),<br/> [simK4_coords.txt](inst/extdata/simK4_coords.txt) (optional) |  [DAPC_simK4_withCoords_v2.pdf](inst/extdata/DAPC_simK4_withCoords_v6.pdf) <br/>[DAPC_simK4_withCoords_v2_BiPlots.pdf](inst/extdata/DAPC_simK4_withCoords_v6_BiPlots.pdf) <br/>[DAPC_simK4_withCoords_v2_densityPlots_PC.pdf](inst/extdata/DAPC_simK4_withCoords_v6_densityPlots_PC.pdf) <br/>[DAPC_simK4_withCoords_v2_densityPlots_DF.pdf](inst/extdata/DAPC_simK4_withCoords_v6_densityPlots_DF.pdf)
```run_sNMF``` | Pipeline for running sNMF (LEA) and graphing results | [simK4.vcf.gz](inst/extdata/simK4.vcf.gz),<br/> [simK4_coords.txt](inst/extdata/simK4_coords.txt) (optional)| [sNMF_simK4_withCoords.pdf](inst/extdata/sNMF_simK4_withCoords.pdf)
```run_fastStructure```  | Pipeline for running fastStructure and graphing results| [simK4.vcf.gz](inst/extdata/simK4.vcf.gz),<br/> [simK4_coords.txt](inst/extdata/simK4_coords.txt) (optional) | coming soon
```runtess``` | Pipeline for running tess3r and generating graphs of results comparable to the other pop structure methods | [simK4.vcf.gz](inst/extdata/simK4.vcf.gz),<br/> [simK4_coords.txt](inst/extdata/simK4_coords.txt)| coming soon
```runeems_snps_setup``` | Generates some of the input files and arranges all necessary inputs in a nice environment for EEMS | coming soon | coming soon



 <!-- - ```eemsgg2raster```-->

### Install ```misc.wrappers``` from GitHub using ```devtools```, ```remotes```, or ```BiocManager```
```
# with devtools:
library(devtools)
devtools::install_github("JeffWeinell/misc.wrappers")

# with remotes:
library(remotes)
remotes::install_github("JeffWeinell/misc.wrappers")

# with BiocManager:
library(BiocManager)
BiocManager::install("JeffWeinell/misc.wrappers")
```

### Software Requirements
Some functions depend on other programs to be installed. In particular, ```run_fastStructure``` requires [fastStructure](https://rajanil.github.io/fastStructure/) to be installed, and ```runeems_snps_setup``` requires [EEMS](https://github.com/dipetkov/eems) to be installed.

#### Download, unpack, and install EEMS
Follow the instructions [here](https://github.com/dipetkov/eems) on how to download the EEMS repository, which includes several programs, including ```runeems_snps``` which you must install from source.

#### Download, unpack, and install fastStructure
Follow the instructions [here](https://rajanil.github.io/fastStructure/) on how to install fastStructure


#### (Optionally) link misc.wrappers to the programs that it depends on. If this is skipped you will need either add the location of these programs to your system path or provide the full paths as arguments when running functions.

If runeems_snps is installed on your PATH, then you can skip this step. Otherwise, use the function ```config_miscwrappers``` to tell misc.wrappers where to find runeems_snps. You only need to do this once.
```
library(misc.wrappers)
# Link the EEMS program 'runeems_snps'
config_miscwrappers(exe.paths="*/eems-master/runeems_snps/src/runeems_snps")

# Link python 2 (required by fastStructure) and the fastStructure program 'structure.py'.
config_miscwrappers(exe.paths=c("*/PATH/TO/python","*/fastStructure-master/structure.py"))

```
### ```run_DAPC``` Run DAPC from SNP data in a VCF file and plot results
```
library(misc.wrappers)
## Example 1: 
# Define path to a VCF file that contains simulated SNP data, or, set the path to a file your own data.
vcf.path    <- file.path(system.file("extdata", package = "misc.wrappers"),"simK4.vcf.gz")
# Run DAPC analyses and save the results to files in your current directory.
run_DAPC(x=vcf.path,format="VCF",kmax=10,samplenames=NULL,reps=30,probs.out=NULL,save.as="DAPC_simK4_withCoords_v4.pdf",include.out=c(".pdf"))

## Example 2: Same SNP dataset as example 1, but here we also provide Lon/Lat coordinates of individuals to geographically interpolate population assignments.
# Path to VCF file
vcf.path    <- file.path(system.file("extdata", package = "misc.wrappers"),"simK4.vcf.gz")
# Path to coordinates file
coords.path <- file.path(system.file("extdata", package = "misc.wrappers"),"simK4_coords.txt")
# Running DAPC analyses
run_DAPC(x=vcf.path, format="VCF", kmax=10, coords=coords.path, samplenames=NULL, reps=30,probs.out=NULL, save.as="DAPC_simK4_withCoords_v4.pdf", include.out=c(".pdf"))
```

### ```run_sNMF``` Run sNMF/LEA from SNP data in a VCF file and plot results
```
library(misc.wrappers)

## Example 1:
# Path to VCF with SNPs
vcf.path    <- file.path(system.file("extdata", package = "misc.wrappers"),"simK4.vcf.gz")
run_sNMF(x=vcf.path,format="VCF",kmax=10,reps=30,save.as="sNMF_simK4.pdf")

## Example 2: Same SNP dataset as example 1, but here we also provide Lon/Lat coordinates of individuals to geographically interpolate admixture coefficients.
vcf.path    <- file.path(system.file("extdata", package = "misc.wrappers"),"simK4.vcf.gz")
coords.path <- file.path(system.file("extdata", package = "misc.wrappers"),"simK4_coords.txt")
run_sNMF(x=vcf.path, format="VCF", coords=coords.path, samplenames=NULL, kmax=10, reps=30, entropy=TRUE, project="new", iter=500, save.as="sNMF_simK4_withCoords.pdf")
```

### ```run_fastStructure``` Run fastStructure from SNP data in a VCF file and and plot results.
```
library(misc.wrappers)
## Example 1:
# Path to VCF with SNPs
vcf.path    <- file.path(system.file("extdata", package = "misc.wrappers"),"simK4.vcf.gz")
# Run fastStructure 30 times each for K=1-10
run_fastStructure(x=vcf.path,kmax=10,reps=30,save.as="fs_simK4.pdf",include.out=c(".pdf"))

## Example 2: Same SNP dataset as example 1, but here we also provide Lon/Lat coordinates of individuals to geographically interpolate admixture coefficients.
vcf.path    <- file.path(system.file("extdata", package = "misc.wrappers"), "simK4.vcf.gz")
coords.path <- file.path(system.file("extdata", package = "misc.wrappers"), "simK4_coords.txt")
run_fastStructure(x=vcf.path, coords=coords.path, kmax=10, reps=30, save.as="fs_simK4_withCoords.pdf", include.out=c(".pdf"))
```

### ```runtess``` Run tess3r from SNP data in a VCF file and plot results
```
library(misc.wrappers)

## Example 1:

```

### ```plot_model``` Create a graphical representation of a demographic model that is defined by a template (.tpl) and estimation (.est) file.
```
library(misc.wrappers)
model.obj  <- plot_model(tpl.path="/PATH/TO/model.tpl",  est.path="/PATH/TO/model.est")
pdf("cartoon_model.pdf",width=10,height=6)
model.obj
dev.off()
```

### ```runeems_snps_setup``` Generate input files for runeems_snps
This function takes as input two files: (1) SNPs in a VCF file and (2) a two-column text file with longitude and latitude coordinates for each individual in the VCF file. Output of this function includes all of the files and directories necessary for running EEMS. Bash files for running EEMS are also generated.

```
library(misc.wrappers)
## Set up input files for runeems_snps
runeems_snps_setup(data="Path/To/SNP/file.vcf",coord="Path/To/LonLat/of/Individuals/file.txt",output.dirpath="Path/To/Directory/That/Doesnt/Exist")
```
Running runeems_snps_setup will create a directory with path output.dirpath and all of the following files:
  - `habitat_outer.pdf`
  - `data/data.diffs`
  - `data/data.outer`
  - `data/data.coords`
  - `mcmc/chain*/`; one for each chain
  - `params/params-chain*.ini`; one for each chain
  - `runeems_snps_chain*.sh`; bash script for running EEMS to generate MCMC chain*

<!---
# I don't remember if this works yet.
### Visualizing results
The function `make_eems_plots` from the reemsplots2 will plot the results, but the maps produced can be difficult to work with further because they are in ggplot objects. The function ``gg2raster`` will convert these to a raster brick object, and optionally save the raster as a geoTIFF file that can be read into GIS software such as QGIS.
```
# Generate the list of ggplot objects
gg <- reemsplots2::make_eems_plots(mcmcpath = `/mcmc/chain*/`)

# Create a raster brick object for each map and save each as a geoTIFF file.
mrates1.brick <- eemsgg2raster(gg.obj=gg$mrates01,file.out="mrates1.tif")
mrates2.brick <- eemsgg2raster(gg.obj=gg$mrates02,file.out="mrates2.tif")
qrates1.brick <- eemsgg2raster(gg.obj=gg$qrates01,file.out="qrates1.tif")
qrates2.brick <- eemsgg2raster(gg.obj=gg$qrates02,file.out="qrates2.tif")
```
--->










