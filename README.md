# misc.wrappers
 **Under Developement**

### Package Description
This package includes functions to fascilitate running some popular programs used for assessing population structure, demographic modelling, and genetic/genomic and phylogenetic/genomic programs, and to summarize results in useful ways. The main functions are shown in the table below along with example input and output files. The values used for function arguments, not shown in table, also control analyses and should be set mindfully - see examples (will add soon) of function usage.

function | example input files | example output files
---|---|---
```plot_model``` | [model10_K3.tpl](inst/extdata/example.tpl), [model10_K3.est](inst/extdata/example.est) | [model10_K3.pdf](inst/extdata/example_model.pdf)
```create.outer``` | [coords.txt](inst/extdata/createouter_exampleInput_coords.txt) | 'method'=1: [outer.txt](inst/extdata/createouter_exampleOutput_method1_outer.text), [outer.pdf](inst/extdata/createouter_exampleOutput_method1_outer.pdf)<br /> 'method'=2: [outer.txt](inst/extdata/createouter_exampleOutput_method2_outer.text), [outer.pdf](inst/extdata/createouter_exampleOutput_method2_outer.pdf)<br /> 'method'=3: [outer.txt](inst/extdata/createouter_exampleOutput_method3_outer.text), [outer.pdf](inst/extdata/createouter_exampleOutput_method3_outer.pdf)
```run_DAPC``` | [simulated_K4.vcf.gz](inst/extdata/simulated_K4.vcf.gz),<br/> *coords.txt* (optional) | 
```run_sNMF``` | [simulated_K4.vcf.gz](inst/extdata/simulated_K4.vcf.gz),<br/> *coords.txt* (optional)| 'coords'=NULL: [simK4_coordsNULL.pdf](inst/extdata/simK4_coordsNULL.pdf)<br/>'coords'="*/coords.txt": [coming soon]
```run_fastStructure```  | [simulated_K4.vcf.gz](inst/extdata/simulated_K4.vcf.gz),<br/> *coords.txt* (optional) | 'coords'=NULL:[simulated_K4_kmax10_fs.pdf](inst/extdata/simulated_K4_kmax10_fs.pdf)<br/>'coords'="*/coords.txt": [coming soon]
```runtess``` | [simulated_K4.vcf.gz](inst/extdata/simulated_K4.vcf.gz),<br/> *coords.txt*| 
```runeems_snps_setup``` | [simulated_K4.vcf.gz](inst/extdata/simulated_K4.vcf.gz),<br/> *coords.txt* | 


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
# Define path to a VCF file that contains simulated SNP data, or, set the path to a file your own data.
vcf_path <- file.path(system.file("extdata", package = "misc.wrappers"),"simulated_K4.vcf.gz")
# Run DAPC analyses and save the results to files in your current directory.
simDAPC <- run_DAPC(vcf="/PATH/TO/SNPs/VCF.vcf",coords="/PATH/TO/LonLat.txt",out="/PATH/FOR/RESULTS.pdf")
# Mostly the same as before, but include a file with coordinates of the samples in the VCF.
# example will be added soon # 
```

### ```run_sNMF``` Run sNMF/LEA from SNP data in a VCF file and plot results
```
library(misc.wrappers)
run_sNMF(vcf="/PATH/TO/SNPs/VCF.vcf",coords="/PATH/TO/LonLat.txt",out="/PATH/FOR/RESULTS.pdf")
```

### ```runtess``` Run tess3r from SNP data in a VCF file and plot results
```
library(misc.wrappers)
runtess(vcf="/PATH/TO/SNPs/VCF.vcf",coords="/PATH/TO/LonLat.txt",out="/PATH/FOR/RESULTS.pdf")
```

### ```run_fastStructure``` Run fastStructure from SNP data in a VCF file and and plot results.
```
library(misc.wrappers)
run_fastStructure(vcf="/PATH/TO/SNPs/VCF.vcf",coords="/PATH/TO/LonLat.txt",out="/PATH/FOR/RESULTS.pdf")
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










