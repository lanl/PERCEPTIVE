# PERCEPTIVE
PERCEPTIVE is an R Shiny based graphical interface for the prediction and annotation of chromatin modifiying enzymes in novel or poorly studied species. PERCEPTIVE applies heuristic models and model species information to make predictions and provides end users with a graphical interpretation of the data.
PERCEPTIVE is broken into two graphical R functions, Pipeline and App. 
The Pipeline function encompasses the entirety of the annotation pipeline (from unassembled genome to predictions). This function will ONLY run in a linux environment, and requires several additional dependencies be downloaded from gdrive.
The App function encompasses the entirety of the data visualization and prediction tool. This tool is agnostic of operating system, and should run on any OS with support for R, and basic CRAN packages.

## Installation
First determine if you are a pipeline user or an app user. Installation will be broken down based on user case.
## Dependencies commond to both functions
Follow these directions irrespective of which function you intend to use.
### Rstudio and R
Both the Pipeline and App functions require a current version of R to be installed, and the R devtools package to download and install PERCEPTIVE from github. 
Here we recommend that users download and install Rstudio (which for some operating systems includes R), found [here](https://posit.co/download/rstudio-desktop/) for easier graphical interfacing, and R version 4.4 (Pile of Leaves). R version 4.4.x can be found for Windows OS [here](https://cran.r-project.org/bin/windows/base/), and Linux and OS X (Darwin) [here](https://cran.r-project.org/src/base/R-4/).
### R devtools
Once you have installed R, start an R terminal using the R GUI, or open RStudio (the terminal can be found at the bottom left corner). Run the following command:
'install.packages("devtools")'.
## Dependencies for just the pipeline function (LINUX only)
[braker.sif](https://drive.google.com/file/d/152hLaqatgFi6k7oyWFv47gTMb_26Sh_j/view?usp=drive_link)
[Perceptivev0.1.sif](https://drive.google.com/file/d/1-44qtlKWFssNO9utKUikWy10yjTFRH7n/view?usp=drive_link)
[Eukaryota.fa](https://drive.google.com/file/d/1WoalwL3oIZfgH7mYAfF0HEbwIhVvFGMM/view?usp=drive_link)

## COPYRIGHT
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
