# misc.wrappers
 Collection of wrappers for programs that I did not create.

### Package Description
This package includes functions to fascilitate running the program runeems_snps. The main function is ```runeems_snps_setup```, which takes as input two files: (1) a SNPs dataset in (Variant-Call Format) VCF and (2) a two-column text file with longitude and latitude coordinates for each individual in the VCF file. Output of this function include all of the files and directories necessary for running EEMS. Bash files for running EEMS are also generated.

### Download, unpack, and install EEMS
Follow the instructions [here](https://github.com/dipetkov/eems) on how to download the EEMS repository, which includes several programs, including ```runeems_snps``` which you must install from source.

### Install ```misc.wrappers``` from GitHub using ```devtools```, ```remotes```, or ```BiocManager```.
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

### Link misc.wrappers to runeems_snps executable
If runeems_snps is installed on your PATH, then you can skip this step. Otherwise, use the function ```config_miscwrappers``` to tell misc.wrappers where to find runeems_snps. You only need to do this once.
```
library(misc.wrappers)
# Provide the full path to the runeems_snps executable
config_miscwrappers(exe.paths="*/eems-master/runeems_snps/src/runeems_snps")
```

### Generate runeems_snps input files
```
library(misc.wrappers)
## runeems_snps_setup function
runeems_snps_setup(output.dirpath="Path/To/Directory/That/Doesnt/Exist",data="Path/To/SNP/file.vcf",coord="Path/To/LonLat/of/Individuals/file.txt")
```
Running runeems_snps_setup will create the output directory defined by output.dirpath, which will contain:
  - `habitat_outer.pdf`
  - `data/data.diffs`
  - `data/data.outer`
  - `data/data.coords`
  - `mcmc/chain*/`; one for each chain
  - `params/params-chain*.ini`; one for each chain
  - `runeems_snps_chain*.sh`; bash script for running EEMS to generate MCMC chain*

### Visualizing results
The function `make_eems_plots` from the reemsplots2 will plot the results, but the maps produced can be difficult to work with further because they are in ggplot objects. The function ``gg2raster`` will convert these to a raster brick object, and optionally save the raster as a geoTIFF file that can be read into GIS software such as QGIS.
```
# Generate the list of ggplot objects
gg <- reemsplots2::make_eems_plots(mcmcpath = `/mcmc/chain*/`)

# Create a raster brick object for each map and save each raster brick as a geoTIFF file.
mrates1.brick <- eemsgg2raster(gg.obj=gg$mrates01,file.out="mrates1.tif")
mrates2.brick <- eemsgg2raster(gg.obj=gg$mrates02,file.out="mrates2.tif")
qrates1.brick <- eemsgg2raster(gg.obj=gg$qrates01,file.out="qrates1.tif")
qrates2.brick <- eemsgg2raster(gg.obj=gg$qrates02,file.out="qrates2.tif")
```











