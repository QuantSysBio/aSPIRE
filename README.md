# aSPIRE
aSPIRE is a tool for the estimation of peptide abundances using label-free quantification mass spectrometry
<img src="aSPIRE_white.png" width="400">

## overview and requirements
*aSPIRE* processes peptide-spectrum matches (PSMs) that were assigned by *inSPIRE*, quantifies them using [Skyline](https://skyline.ms/project/home/begin.view) and constructs a generation kinetic for each identified peptide.
It does not do any heavy computing, thus can be run on a **Linux, Mac or Windows laptop** on a single CPU. (**Note**, that software testing was done only on macOS/Linux).

The main steps of *aSPIRE* are:
- parsing of *inSPIRE* assignments (PSMs), peptide mapping and creation of input files for Skyline
- quantification of peptides across all raw files provided using Skyline (can be run either automatically by using a [Docker image](https://hub.docker.com/r/chambm/pwiz-skyline-i-agree-to-the-vendor-licenses) or manually on a Windows machine, see below for more details)
- aggregation of peptide abundances to generation kinetics
- removal of potential synthesis errors and contaminants
- data normalisation and filtering of noisy data points
- results visualisation (peptide generation kinetics, total ion chromatograms, peptide coverage maps, residue maps etc)

## installation of dependencies
### 
### Skyline


## user input

## execution

## output

