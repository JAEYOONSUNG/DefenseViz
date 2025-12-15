# DefenseFinderViz
![DefenseFinder_Heatmap](https://github.com/user-attachments/assets/68e953c1-d568-4b56-89f2-1107a60a6f6c)

## Project Introduction

Welcome to our GitHub repository, where we're excited to share a series of workflows designed to streamline processes in systems biology. This repository is composed of various scripts, each tailored to specific tasks within our broader research framework. Additionally, we're providing access to a curated database to enhance your research capabilities.
- Principal investigator: Dong-Woo Lee
- Project lead: Jae-Yoon Sung
- Maintainers: Jae-Yoon Sung
- Contributors: Jae-Yoon Sung

## Getting Started:
To begin using our resources, please follow the steps outlined in our documentation. 
Whether you're looking to integrate our scripts into your existing projects or explore our database for new insights, we've provided all the necessary instructions to get you started.

## Installation
### Requirements

The DNMB is supported for macOS, Linux and Windows machines, which can provide an environment for using R.
It requires R version >=4.2.1 for release, and R version >=4.3 for devel.

One of the third-party functionalities is not available for Windows and MacOS machines (InterProScan).

The [EggNOG-mapper webserver](http://eggnog-mapper.embl.de), allows users to input sequences in FASTA format based on locus_tag identifiers and receive results in either XLSX or CSV format. Additionally, the standalone version available on GitHub is compatible with DNMB.

InterProScan requires a Linux operating system. Without access to Linux, you can proceed with the analysis up to Eggnog-mapper in the annotation stage, but you won't be able to obtain information about motif analysis.


To download and install R, see the [R-project website](https://www.r-project.org/).

To download and install InterProScan, see the [InterProScan github](https://github.com/ebi-pf-team/interproscan).

#### Warning
The basic file for genomic analysis, known as a GenBank file, requires both sequence and annotation in full-format files such as gbff, gb, or gbk. Additionally, GenBank prefers a format based on the GeneMarkS2+ pipeline, and using a different annotation pipeline to obtain GenBank files may lead to errors.

## Anaylsis flow
## Prerequisites
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "ComplexHeatmap", "Peptides"))

install.packages(c("qdap", "seqinr", "stringr", "stringi", "splitstackshape", "gtools", "ggplot2", "ggseqlogo", "circlize", "grid", "gridExtra","plyr", "dplyr", "tidyr", "readr", "reshape2", "data.table", "tibble", "qdap", "openxlsx"))
```
- **Note:** If you encounter issues installing the qdap package, try installing it with the following command:
```r
install.packages("qdap", INSTALL_opts = "--no-multiarch")
```


- **Note:** Java and the rJava package must be installed and configured to enable .xlsx output using this package.
1. Install Java Development Kit (JDK):
        Download and install the appropriate JDK for your operating system from the Oracle website or OpenJDK.
2. Install rJava Package in R:
```r
install.packages("rJava")
library(rJava)
```

3. Set $JAVA_HOME Path:

You need to set the environment variable JAVA_HOME to point to the location of your JDK installation.
	
 •	On Windows:

1.	Install Java jdk (https://www.oracle.com/kr/java/technologies/downloads/)
2.	Check the “System Variables,” :
	•	Variable name: JAVA_HOME
	•	Variable value: The path to your JDK installation 

```bash
# Print the current value of the JAVA_HOME environment variable.
echo %JAVA_HOME%  #(e.g., C:\Program Files\Java\jdk-18)

# Set the JAVA_HOME environment variable to point to the Java Development Kit (JDK) installation.
# Replace [version] with your installed JDK version (e.g., jdk-18).
setx JAVA_HOME "C:\Program Files\Java\jdk[version]"

# Update the system PATH to include the bin directory of the JDK.
setx PATH "%JAVA_HOME%\bin;%path%"

# Check the installed Java version to confirm that the correct version is being used.
JAVA -version
```
3.	Restart R or RStudio.
	
 •	On macOS/Linux:
-If Xcode is not installed, you may encounter compiler issues during package installation. To resolve this, install Xcode from the App Store.
Add this line to your .bash_profile or .bashrc (depending on the shell):
```bash
# Navigate to your Java installation directory to check available Java versions
/Library/Java/JavaVirtualMachines/[my_java_folder]/Contents/Home # check my java list

# Open your .bash_profile (or .bashrc) file for editing
vi ~/.bash_profile # edit bash profile

# Press 'i' to enter insert mode in the vi editor
i # insert mode

# Add or update the JAVA_HOME environment variable with the path to your Java installation
export JAVA_HOME=/Library/Java/JavaVirtualMachines/[my_java_folder]/Contents/Home
# Add Java's bin directory to the system PATH variable so that Java commands can be run from the terminal
export PATH=${PATH}:$JAVA_HOME/bin

# Save the changes and exit the vi editor. ":wq!" means "write" (save) and "quit" (exit) forcefully
: # activate command line
wq! # save

# Apply the changes made to the .bash_profile or .bashrc immediately (without needing to restart the terminal)
source ~/.bash_profile  ## or ~/.bashrc #apply changes

# Verify that JAVA_HOME is set correctly by printing its value
echo $JAVA_HOME # validation
```

## Install DNMB R package
```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("JAEYOONSUNG/DNMB")
```

# Quick start
## Run DNMB analysis
```r
setwd([GenBank directory]) # Set the working directory to the location where your GenBank files are stored.
library(DefenseFinderViz)
```
---
**DefenseFinderViz**


```r
DefenseFinder_Heatmap()
```
- **Note:** [Strain of interest].defense_finder_systems.tsv output are used for merging data. GenBank’s SOURCE field is used for extracting names.
- **Note:** protien coding sequence (.faa) output was used for defense-finder analysis (https://github.com/mdmparis/defense-finder)

## Contributing
We welcome contributions from the community! If you have suggestions for improvements, additional scripts, or updates to the database, please see our contributing guidelines for more information on how to get involved.


## License
This project is released under MIT licence, which allows for both personal and commercial use, modification, and distribution of our work, provided that proper credit is given.

We hope our resources will prove invaluable to your research in systems biology. For any questions or feedback, please don't hesitate to reach out through our GitHub issues or contact section.

## Citation
If you use this piepline, please cite:
```
[DNMB] DNMB: A Strategic Blueprint for the Domestication of Geobacillus stearothermophilus as a Thermophilic Platform using the DNMB Suite.
             Jae-Yoon Sung, Mun Hoe Lee, Hyungbin Kim, Dariimaa Ganbat, Hyun-Woo Cho, Sang Jae Lee, Seong Bo Kim, and Dong-Woo Lee. 2024.
             XXX, XXX, https://doi.org/XXX
```
Please, cite also the underlying algorithm if it was used for the search step of DNMB:
```
[DefenseFinder] DefenseFinder: Systematic and quantitative view of the antiviral arsenal of prokaryotes. Tesson, F., Hervé, A., Mordret, E. et al. 
				Nat Commun 13, 2561 (2022). https://doi.org/10.1038/s41467-022-30269-9 
```
