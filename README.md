# misc.wrappers
 Collection of wrappers for programs that I did not create.

### Package Description
This package includes functions to fascilitate running the program runeems_snps. The main function is ```runeems_snps_setup```, which takes as input two files: (1) a SNPs dataset in VCF-format and (2) a two-column text file with longitude and latitude coordinates for each individual in the VCF file. Output of this function include all of the files and directories necessary for running EEMS. Bash files for running EEMS are also generated.

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

### Setup misc.wrappers 
If runeems is installed on your system PATH, then you can skip this step. Otherwise, use the function ```config_miscwrappers``` to permanently tell misc.wrappers where to find EEMS executables.
```
library(misc.wrappers)


```

