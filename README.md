# MeristemBasic: 3D dynamics of Arabidopsis thaliana WUS/CLV3 network

MeristemBasic is a hybrid ODE/stochastic based signaling simulator that tracks the changes in concentration of different effectors of the WUS/CLV3 network over time under user defined conditions. Several utilties are provided to prepare and process data both before and after simulation. 

## Getting Started

### Dependencies and requirements

MeristemBasic has been tested to run natively on Ubuntu Linux and Windows through Windows Subsystem for Linux. 


### Installing

Linux

unzip repository to directory of your choice


Windows

Install Windows Subsystem for Linux
https://docs.microsoft.com/en-us/windows/wsl/install-win10 (You just need WSL 1 although WSL 2 is also good) 
On updated systems WSL can be installed with a single command through Powershell
wsl --install
Restart system to complete installation
Open Powershell or open and run Ubuntu through the Start Menu
Download the repository from Bitbucket and extract it.



### Sample Usage

**Linux**

Requirements: current versions of python and R as well as libraries listed in the scripts

Open Terminal and browse to repo debug build directory

sample run command
./MeristemBasic_p parameterlist062021FcAcWt.model StandardSteadyStateWUSGradient01272A.init rk6 > results062021FcAcWt

results files will be generated in the directory (in this case result062021FcAcWt and extraot062021FcAcWt)

#Preparing the Data

place results files in main repo folder

take the python script MerSimProcessor and mersimconfig from the scripts folder and place them in the main repo directory. Update mersimconfig so that the file paths match your system.

Run the script (eg python MerSimProcessor). Pick option 6 

The partitioned results should be output to the Results/Raw_DataFiles folder

To visualize the results take the last file from the output results folder and rename it EWT.csv or run a terminal command such as 
cp ./Results/Raw_DataFiles/extraresults062021FcAcWt/extraresults099.csv ./EWT.csv 

to copy it to the main directory. 

#Visualizing the Data

Start R

Ensure that your working directory in R is set to where the datafiles are located

open and run the scripts (in the script folder)cvsShellerGeneric2.R (change the working directory in this script) and MerBasicGrapher.R 


**Windows Subsystem for Linux (Windows 10)**


Open Powershell and navigate to the extracted repository 

change to wherever the directory is on your system) using the cd command

run the Simulator as you would a local Linux install but with wsl at the beginning if you are using Powershell
example

wsl ./MeristemBasic_p parameterlist062021FcAcWt.model StandardSteadyStateWUSGradient01272A.init rk6 > results062021FcAcWt



## Acknowledgement
The base simulator code is forked from Organism Simulator by Henrik Jonsson used in [1]

[1] Yadav RK, Perales M, Gruel J, et al. Plant stem cell maintenance involves direct transcriptional repression of differentiation program. Mol Syst Biol. 2013;9:654. doi:10.1038/msb.2013.8
