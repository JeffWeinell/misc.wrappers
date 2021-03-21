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

The output of ```runeems_snps_setup``` includes the following:
  - `output.dirpath/`
    - `habitat_outer.pdf`
    - `data/`
      - `data.diffs`
      - `data.outer`
      - `data.coords`
    - `mcmc/`
      - `XXXDemes-chain1/`
      - `XXXDemes-chain*/`
    - `params/`
      - `params-chain1_XXX.ini`
      - `params-chain*_XXX.ini`




  - `*.ini` files: parameter settings for each chain
  - `data.diffs`: file containing the pairwise matrix of differences (distances)
  - `data.outer`: file containing the coordinates defining the habitat region in which to estimate effective migration.
  - `data.coords`: a copy of the input coordinates defining the location of samples in the VCF
  `habitat_outer.pdf`: a pdf image showing the the region 
  - 

### Example 










