# PERCEPTIVE
PERCEPTIVE is an R Shiny based graphical interface for the prediction and annotation of chromatin modifying enzymes in novel or poorly studied species. PERCEPTIVE applies heuristic models and model species information to make predictions and provides end users with a graphical interpretation of the data.
PERCEPTIVE is broken into two graphical R Shiny functions: the PERCEPTIVE app function and the Pipeline function. 

The PERCEPTIVE app function encompasses the entirety of the data visualization and prediction tool. This tool is agnostic of operating system, and should run on any OS with support for R and basic CRAN packages.

The Pipeline function encompasses the entirety of the annotation pipeline (from unassembled genome to predictions). This function will ONLY run in a Linux environment, and requires several additional dependencies be downloaded from gdrive.

The manuscript preprint which details this package can be found [here](https://www.biorxiv.org/content/10.1101/2024.11.20.624555v1)

## Installation
First determine if you are an app user (you want to view data which has already been processed by another user) or a pipeline user (you intend to process raw data for subsequent visualization). Installation is broken down based on user case below.
## Dependencies common to both functions
Follow these directions irrespective of which function you intend to use.
### Rstudio and R
Both the Pipeline and App functions require a current version of R to be installed, and the R devtools package to download and install PERCEPTIVE from github. 
Here we recommend that users download and install Rstudio (which for some operating systems includes R), found [here](https://posit.co/download/rstudio-desktop/) for easier graphical interfacing, and R version 4.4 (Pile of Leaves). R version 4.4.x can be found for Windows OS [here](https://cran.r-project.org/bin/windows/base/), and Linux and OS X (Darwin) [here](https://cran.r-project.org/src/base/R-4/).
### R devtools
Once you have installed R, start an R console using the R GUI or open RStudio (the console can be found at the bottom left corner). Run the following command in the R console:
```
install.packages("devtools")
```
### To install the PERCEPTIVE package
Run the following command in the R console:
```
devtools::install_github("https://github.com/lanl/PERCEPTIVE.git", dependencies = TRUE)
```
Be certain to accept "All" dependencies if prompted

## Dependencies for just the pipeline function (LINUX only)
The Pipeline function requires Singularity 4.2.x or greater, two Singularity images and one OrthoDB database to run properly.
### Install Singularity (Apptainer)
To install Singularity run the following in a terminal:

Install GO compiler
```
wget https://go.dev/dl/go1.23.3.linux-amd64.tar.gz && \
sudo tar -C /usr/local -xzvf https://go.dev/dl/go1.23.3.linux-amd64.tar.gz && \
rm https://go.dev/dl/go1.23.3.linux-amd64.tar.gz
```
Setup GO Environment
```
echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc
```
Install Singularity
```
mkdir -p $GOPATH/src/github.com/sylabs && \
cd $GOPATH/src/github.com/sylabs && \
wget https://github.com/sylabs/singularity/releases/download/v4.2.1/singularity-ce-4.2.1.tar.gz && \
tar -xzf singularity-${VERSION}.tar.gz && \
cd ./singularity && \
./mconfig && \
make -C ./builddir && \
sudo make -C ./builddir install

```
Alternatively, visit the Singularity documentation [here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) for more information. Note, as of December 2024 the Singularity documentation is out of date, and provides instructions to install Singularity version 3.x. PERCEPTIVE requires 4.2.x or greater.

### Download Singularity images and OrthoDB database.
PERCEPTIVE requires several singularity images and a pre-compiled database for BRAKER3. These are provided via gdrive at the links below. On startup PERCEPTIVE will prompt the user for paths to these files. 

The BRAKER3 image can be generated following the instruction in the [BRAKER3 documentation](https://github.com/Gaius-Augustus/BRAKER), or can be downloaded from gdrive [here](https://drive.google.com/file/d/152hLaqatgFi6k7oyWFv47gTMb_26Sh_j/view?usp=drive_link).

The Perceptive Singularity image can be downloaded from gdrive [here](https://drive.google.com/file/d/1-44qtlKWFssNO9utKUikWy10yjTFRH7n/view?usp=drive_link).

The pre-partitioned OrthoDB database compiled for BRAKER3 can be downloaded from gdrive [here](https://drive.google.com/file/d/1WoalwL3oIZfgH7mYAfF0HEbwIhVvFGMM/view?usp=drive_link). Newer [versions](https://github.com/Gaius-Augustus/BRAKER) may be available from the maintainers of BRAKER3, but support is not guaranteed.

PERCEPTIVE also requires a license file for GeneMark, which can be found [here](https://genemark.bme.gatech.edu/license_download.cgi). Please select: GeneMark-ES/ET/EP+ ver 4.72_lic,  LINUX 64 kernel 3.10 - 5. Upon download the genemark license will be gzipped. PERCEPTIVE expects the license to be unzipped. 
To unzip run the following:
```
gzip -d [pathtoyourlicensefile]
```
## Usage
As above usage will be broken down by for the two functions PERCEPTIVE and Pipeline. For the time being the graphical interface must be started from the command-line. A launcher application may become available in the future.

### PERCEPTIVE Usage
To run the PERCEPTIVE Graphical Interface for interpretation and prediction of chromatin modifying enzymes open an R terminal and run the following:
```
PERCEPTIVEv3::PERCEPTIVE()
```
The graphical interface will prompt you for files generated by the Pipeline. These will be found within the sub directory FilesForGUI (simply select the directory), which can be found within the species directory generated by the Pipeline function, at the location specified when the Pipeline was run.

### Pipeline Usage
To run the Pipeline Graphical Interface for annotation of chromatin modifying enzymes open an R terminal and run the following:
```
PERCEPTIVEv3::Pipeline()
```
### Options
The Pipeline function will prompt the user for several inputs prior to running. These input represent four run modes (although not explicitly referred to as modes 1-4 within the interface). The first two modes are recommended, and in most cases mode one is sufficient. Modes three and four are experimental. Modes three and four assume a "typical" diploid genome. In most cases it would be preferable to assemble a de novo genome before running through the Pipeline.
#### One: Assembled Genome without RNA-seq
In this mode the input is solely an unmasked assembled genome. This genome could be generated in-house, or downloaded from JGI, NCBI, ENSEMBLE, etc. In this mode BRAKER3 will use the Eukaryotic OrthoDB database for protein prediction. (Selections: Would you like to perform de novo genome assembly: no, RNA-seq: no)
#### Two: Assembled Geneome with RNA-seq
In this mode the input is an unmasked assembled genome and unaligned RNA-seq reads from a total RNA library preparation (note: sequencing depth must be at minimum 20X coverage for the transcriptome). This mode may perform better if a species is highly esoteric, ie., it is expected to possess many proteins with high divergence from other Eukaryotes. (Selections: Would you like to perform de novo genome assembly: no, RNA-seq: yes)
#### Three: Genome assembly without RNA-seq
Similar to mode one, this mode solely requires input of either unassembled short or long reads. In this mode de novo genome assembly will be performed, and the assembly generated will be used as input to the Pipeline. The user must also provide an approximate estimated genome size (for Canu). We recommend using jellyfish or GenomeSource for estimation of genome size. Alternatively, if a close relative for your species has been assembled, the size of that genome can be used as an estimate. (Selections: Would you like to perform de novo genome assembly: yes, short(velvet) or long (canu), RNA-seq: no)
#### Four: Genome assembly with RNA-seq
Similar to mode three, this  mode requires input of either unassembled short or long reads. In this mode de novo genome assembly will be performed, and the assembly generated will be used as input to the Pipeline. The user must also provide an approximate estimated genome size (for Canu). We recommend using jellyfish or GenomeSource for estimation of genome size. Alternatively, if a close relative for your species has been assembled, the size of that genome can be used as an estimate. This mode additionally requires unaligned RNA-seq data from a total RNA library preparation, (see note above). (Selections: Would you like to perform de novo genome assembly: yes, short (velvet) or long (canu), RNA-seq: yes)
#### Options three and four further details
*Velvet: If Velvet is to be used with short reads, either single or paired end reads can be used. Please specify as appropriate.

*Canu: If the user so chooses, additional arguments can be passed to Canu within PERCEPTIVE. This has not been exhaustively tested, and the user is warned in advance. Please monitor the terminal output for appropriate error and warning messages from Canu and PERCEPTIVE.


### Aditional Option: Blast to database of user choice (override PERCEPTIVE defaults)
To override the PERCEPTIVE model organism blast database and blast against a user-defined database, follow these directions to generate a database:

Download reference files specfic to organisms of interest, identified by [taxon number](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/summary/taxonomy/datasets_summary_taxonomy_taxon/) using NCBI [datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) tool:
```
datasets download genome taxon [your taxa here] --reference --include protein
```
Unzip and extract relevant .faa files. **If you want to create a database containing multiple organisms, the .faa (protein) files can be concatenated.**
```
unzip ncbi_dataset.zip
cat *.faa > total.faa
```
The concatenated files must now be run through [segmasker](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) from blast+.
```
segmasker -in total.faa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out totalprot.asnb
```
Finally, using makeblastdb from [blast+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) you can make a blast database. The resulting folder will be used by PERCEPTIVE instead of the included database of model organism sequences. It is recommended to retain the naming below for directories.
```
makeblastdb -in total.faa -title modelorgs -mask_data totalprot.asnb -dbtype prot -out modelorgsprot -parse_seqids
mkdir databaseprot/
mv *modelorgsprot* databaseprot/
```

## Errors
A 'common' error occurs when using GeneMark-ES with _poorly_ assembled genomes. In the case where all or most contigs are very short (<50kB) GeneMark will fail to run providing no gene predictions. PERCEPTIVE will attempt to salvage this by lowering the threshold for contig size to 10kB. However, in most cases contigs smaller than this represent a poor quality assembly, and resequencing, or reassembly with Velvet or Canu is suggested. The hallmark of this error is the warning message: GeneMark fails! seen in the terminal. Additional error handling will be added at a future date.

## Test Data
A folder including all files neccessary to test the PERCEPTIVE app function is included [here](https://drive.google.com/drive/folders/1nhPLN9m4bMLwtM_H0rF81UbEbCayKo1t?usp=drive_link) on gdrive for Nannochloropsis gaditana. Download all files and place in a single folder. This folder will be the folder selected within the PERCEPTIVE app GUI.
To test the Pipeline function, it is suggested to run the Pipeline in mode one using the Nannochloropsis gaditana genome for CCMP 1894 [here](https://genome.jgi.doe.gov/portal/pages/accessDenied.jsf?state=%27anonDownload%27) from JGI. At a later date GEO or SRA accessions may be made available with RNA-seq data, long read data, and short read data for testing, but is not currently provided.


## COPYRIGHT
This package relies on several dependencies and pipelines. The authors make no claim of ownership, copyright, or maintenance for any or all of these dependencies.

DOE NNSA O# (O4768) 

© 2024. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

(End of Notice)

This program is Open-Source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

(End of Notice)
