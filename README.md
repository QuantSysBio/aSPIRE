# aSPIRE
aSPIRE is a tool for the estimation of peptide abundances using label-free quantification mass spectrometry
<img src="aSPIRE_white.png" width="400">

## Overview and requirements
*aSPIRE* processes peptide-spectrum matches (PSMs) that were assigned by *inSPIRE*, quantifies them using [Skyline](https://skyline.ms/project/home/software/Skyline/begin.view) and constructs a generation kinetic for each identified peptide.
It does not do any heavy computing, thus can be run on a **Linux, Mac or Windows laptop** on a single CPU. However (**Note**, that software testing was done only on macOS/Linux).

The main steps of *aSPIRE* are:
- parsing of *inSPIRE* assignments (PSMs), peptide mapping and creation of input files for Skyline
- quantification of peptides across all raw files provided using Skyline (can be run either automatically by using a [Docker image](https://hub.docker.com/r/chambm/pwiz-skyline-i-agree-to-the-vendor-licenses) or manually on a Windows machine, see below for more details)
- aggregation of peptide abundances to generation kinetics
- removal of potential synthesis errors and contaminants
- data normalisation and filtering of noisy data points
- results visualisation (peptide generation kinetics, total ion chromatograms, peptide coverage maps, residue maps etc)

## Installation of dependencies
The following instructions need to be **executed oncy once to set up *aSPIRE***. Once this is done, you can directly progress to the execution for any further runs.

### Snakemake and Conda
*aSPIRE* relies on [Conda](https://docs.conda.io/en/latest/) and Snakemake.
In order to install Conda, click on this [link](https://docs.conda.io/en/latest/miniconda.html) and follow the installation guidelines for your respective operating system.  
After installing Conda, you need to install Snakemake. The Snakemake installation procedure is described [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Briefly, open the terminal on your computer and paste the following lines sequentially:

    conda install -n base -c conda-forge conda
    conda activate base
    conda create -c conda-forge -c bioconda -n aspire snakemake

Additionally, you might need to run `conda update conda`. We repeatedly faced some issues around Conda installations, particularly with devices with Apple M1 cores. Please refer to [this tutorial](https://pad.gwdg.de/s/7C3rWC3w2#) for troubleshooting your installation.

Download this repository as a .zip file (click on *Code* at the upper right corner of this repository --> Download ZIP), move the .zip file in the desired directory on your computer and unpack it.
Open the terminal in this directory and enter:

    conda activate aspire

### Skyline
For the extraction of peptide abundances, *aSPIRE* relies on [Skyline](https://skyline.ms/project/home/software/Skyline/begin.view). There are two ways to execute *aSPIRE* and Skyline.
1. If you have limited experience with command line, a Windows machine available and only a few samples to process we recommend **manual execution**.
2. For high-throughput data processing with *aSPIRE* we recommend executing Skyline via **command line** in automated fashion.

Choose one of these methods according to your requirements and follow the respective instructions below.

#### Manual execution
Install Skyline on your Windows machine following the [vendor's instructions](https://skyline.ms/project/home/software/Skyline/begin.view). In the config file, remember to change the flag

    automatedSkyline: no

to `no` (see [config file](#config-file) section for detailed instructions). Make sure to follow the instructions in `Skyline_tutorial.pdf` during execution.

#### Command line execution
The command line options requires the installation of [Docker](https://www.docker.com/). Skyline is containerised in [this](https://hub.docker.com/r/chambm/pwiz-skyline-i-agree-to-the-vendor-licenses) Docker image. Pull the image and make sure that it works:

    docker pull chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
    docker run -it --rm chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine SkylineCmd --help

You should see information for all possible flags printed. In the config file, remember to change the flag

    automatedSkyline: yes

to `yes` (see [config file](#config-file) section for detailed instructions).

## User input
All input information should be provided in the `data/` folder. Please have a look at the files `config.yaml` and `sample_list.csv` in that folder.

### config file
Open `data/config.yaml` using a text editor such as VS code or Sublime. Fill in the information as in the example:

    protein_name: CaM
    spAngle: 0
    qVal: 0.01
    RT: 150000000
    raw_file_loc: /absolute_path/to/folder_with_all_raw_files/
    automatedSkyline: yes
    skyline_report: MS1_HPR

- `protein_name`:
- `spAngle`: Do not change unless you have a reason to do so.
- `qVal`: False discovery rate for which PSMs are filtered. Set to 1% per default.
- `RT`: Do not change unless you have a reason to do so.
- `raw_file_loc`: Absolute path pointing towards a folder that contains all raw files.
- `automatedSkyline`:
- `skyline_report`: do not change

### sample list
Open `data/sample_list.csv` with MS Excel/Numbers etc or a text editor. Fill in the information as in the example:

| protein_name | substrateID | substrateSeq | digestTime | biological_replicate | final_assignments | raw_file |
| ----- |  ----- |  ----- |  ----- |  ----- |  ----- |  ----- |
CaM	| CaM	| CaM.fasta| 0 | 	G1	| IL37b_finalAssignments.csv | WSoh_101121_151121_HFGoe_G1_IL37b_0h_R1.raw
CaM	| CaM	| CaM.fasta	| 4	| G1	| IL37b_finalAssignments.csv	| WSoh_101121_151121_HFGoe_G1_IL37b_4h_R1.raw
CaM	| CaM	| CaM.fasta	| 0	| G2	| IL37b_finalAssignments.csv	| WSoh_101121_151121_HFGoe_G2_IL37b_0h_R1.raw
CaM	| CaM	| CaM.fasta	| 2	| G2	| IL37b_finalAssignments.csv	| WSoh_101121_151121_HFGoe_G2_IL37b_2h_R2.raw

A few general remarks:
- Please always provide a full kinetic including the zero hours / no proteasome control measurements.
- You can put any number of proteins / experiments in the sample list. However, **only a single protein can be processed at once**. Specify the `protein_name you` would like to analyse in the config file (see above).
- Do not put any white spaces, umlauts or empty rows in the sample list. Also, every column has to be filled, even if there are duplicated entries (*e.g.*, substrateID).

**Explanation of sample list columns:**
- `protein_name`

After filling in all information, save the sample list in .csv format. **Important!** In case you are creating/editing the `sample_list.csv` file in Microsoft Excel, make sure to save it as actual comma-separated file. *I.e.*, Save as --> Comma-separated Values (.csv) To check that the sample list is in the correct format, you can open it using a text editor and verify that the columns are separated by commas (and NOT semicolons).

## Execution

## output

